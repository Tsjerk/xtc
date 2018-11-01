#!/usr/bin/env python

# Fast operations on a Gromacs XTC trajectory
# (c)2013 Tsjerk A. Wassenaar

from struct import unpack
import sys
import os


# Very simple option class
class Option:
    def __init__(self,func=str,num=1,default=None,description=""):
        self.func        = func
        self.num         = num
        self.value       = default
        self.description = description
    def __nonzero__(self): 
        if self.func == bool:
            return self.value != None
        return bool(self.value)
    def __str__(self):
        return self.value and str(self.value) or ""
    def setvalue(self,v):
        if len(v) == 1:
            self.value = self.func(v[0])
        else:
            self.value = [ self.func(i) for i in v ]

FNM = []
selection = []

# Description
desc = "Fast operations on Gromacs XTC files."

# Option list
options = [
#   option                  type number default description
    ("-f",        Option(FNM.append,  1,        None, "Input XTC file(s)")),
    ("-o",        Option(str,         1,   "out.xtc", "Output XTC file")),
    ("-reverse",  Option(bool,        0,        None, "Reverse trajectory")),                               #
    ("-last",     Option(bool,        0,        None, "Extract last frame(s)")),                            #
    ("-cat",      Option(bool,        0,        None, "Concatenate trajectories, keeping double frames")),  #
    ("-stitch",   Option(bool,        0,        None, "Concatenate trajectories, removing double frames")), #
    ("-select",   Option(selection.append,  1,  None, "Select frames")),                                    #
    ("-info",     Option(bool,        0,        None, "Print trajectory information")),                     #
    ("-times",    Option(bool,        0,        None, "Extract times of all frames")),                      #
    ("-steps",    Option(bool,        0,        None, "Extract steps of all frames")),                      #
    ("-split",    Option(bool,        0,        None, "Split trajectory in individual frames")),            #
#    ("-dump",     Option(str,         1,        None, "Dump the frame closest to the time, frame or step, specified as, e.g.,")),
#"%12s%s"%("","t=50 (time in ns), f=10 (frame number), s=1000 (step number)"),
    ]


args = sys.argv[1:]
if all([i.endswith(".xtc") for i in args]):
    args = ["-info"]+[i for j in args for i in ("-f", j)]


if '-h' in args or '--help' in args:
    print "\n",__file__
    print desc or "\nSomeone ought to write a description for this script...\n"
    for thing in options:
        print type(thing) != str and "%10s  %s"%(thing[0],thing[1].description) or thing
    print
    sys.exit()


# Convert the option list to a dictionary, discarding all comments
options = dict([i for i in options if not type(i) == str])


# Process the command line
while args:
    ar = args.pop(0)
    options[ar].setvalue([args.pop(0) for i in range(options[ar].num)])


###


# Concatenate files. Simplest case (could as well use 'cat')
if options["-cat"] and not options["-reverse"]:
    with open(options["-o"].value,"w") as out:
        for fnm in FNM:
            with open(fnm) as x:
                out.write(x.read())
    sys.exit()


###


def r(a):
    """ Make a slice from a string like start:stop:step (default 0:-1:1). """
    a = a.split(":")
    if len(a) == 1:
        i = int(a[0])
        return slice(i,i+1,1)
    x = a[0] and int(a[0]) or 0
    y = a[1] and int(a[1]) or None
    s = len(a) == 3 and int(a[2]) or 1
    return slice(x, y, s)


def strseek(stream,tag,bufsize=10000):
    """ Find every position in file where tag is found. """
    v = len(tag)
    x = stream.read(bufsize)
    n = 0
    while len(x) >= v:
        m = x.find(tag)
        if m > -1:
            # If we found the tag, the position is the total length 
            # read plus the index at which the tag was found 
            n += m
            yield n
            # Now we strip the part up to and including the tag
            x = x[m+v:]    
            n += v
        elif len(x) > v:
            # If the tag was not found, strip the string, keeping only
            # the last v-1 characters, as these could combine with the
            # next part to be read in.
            # n is incremented with the number of characters taken from
            # the current string x (len(x)-v+1)
            n += len(x)-v+1
            x = x[1-v:]
        if len(x) <= v:
            x += stream.read(bufsize)


class XTC:
    def __init__(self, fnm, last=False):
        """Open the file and read the first frame."""

        self.name  = fnm
        self.file  = open(fnm,"rb")
        header     = self.file.read(92)                      # The XTC header, including box, time, etc.
        self.tag   = header[:8]                              # Tag: magic number and number of atoms

        # XTC:
        #    What                 Bytes Format
        #  0. Magic number (4)    0 -  4 l
        #  1. Atoms (4)           4 -  8 l
        #  2. Step  (4)           8 - 12 l 
        #  3. Time (4)           12 - 16 f
        #  4. Box (9*4)          16 - 52 fffffffff
        # 13. Atoms (4)          52 - 56 l           (checked to be equal to b.)
        # 14. Precision (4)      56 - 60 f
        # 15. Extent (6*4)       60 - 84 llllll      (MIN: x,y,z, MAX: x,y,z)
        # 21. smallidx (4)       84 - 88 l
        # 22. Size in bytes (4)  88 - 92 l 
        # 23. Coordinates (compressed)
        unpacked   = unpack(">lllfffffffffflfllllllll", header)
        self.atoms = unpacked[1]                             # Number of atoms
        self.step1 = unpacked[2]                             # Step of first frame
        self.start = 0.001*unpacked[3]                       # Time of first frame
        self.prec  = unpacked[14]                            # Precision of coordinates
        self.size1 = unpacked[22]+92                         # Size of frame in bytes (including header)
        self.frame = header + self.file.read(self.size1)     # The whole frame

        self.total = self.file.seek(0,2) or self.file.tell() # Total size of the trajectory

        self.selection = None
        self.times     = []
        self.steps     = []
        self.frames    = []
        self.nframes   = 0
        self.lengths   = []

        if last:
            self.lastFrames()


    def __repr__(self):
        stuff = (self.name,self.step1,self.step2,self.start,self.end,self.nframes,self.atoms,self.total)
        return "%-20s: steps %10d to %10d, from %10.5f to %10.5f ns, %5d frames, %5d atoms, total size %d"%stuff


    def __iter__(self):
        if not self.frames:
            self.allFrames(selection=self.selection)

        for i in range(len(self.frames)):
            # Wind to the start of this frame
            self.file.seek(self.frames[i])
            # Return the frame:
            yield self.times[i], self.file.read(self.lengths[i])
            

    def lastFrames(self):
        """ Read the last frame """

        # Jump to somewhere before last frame
        # Do check whether the trajectory consists of one frame...
        try:
            self.file.seek(-3*self.size1/2,2)
        except IOError:
            # Ended up before file; Only one (complete) frame here.
            # Rewinding to start
            self.file.seek(0)

        here       = self.file.tell()
        frame      = self.file.read()
        where      = frame.index(self.tag)
        self.frame = frame[where:]
        self.step2 = unpack(">l",self.frame[8:12])[0]
        self.end   = 0.001*unpack(">f",self.frame[12:16])[0]
        self.size2 = 92 + unpack(">l", self.frame[88:92])[0]

        # The XTC compression may cause variable frame sizes
        # The first and last frames should be close to an 
        # upper and lower bound of the frames size, assuming
        # that the size settles during equilibration.            
        # Using average frames size.
        # There are often some padding bytes;
        # we add one to the length (should check this),
        # and we round down. The number may well be a bit off.
        self.nframes = int(1.0*self.total/(0.5*self.size1+0.5*self.size2+1))


    def allFrames(self, selection=None, steps=False, times=False, reverse=False):
        """ Get all frame positions and sizes. """

        # Rewind file
        self.file.seek(0)
        
        # Find frames
        self.frames  = [ i for i in strseek(self.file,self.tag) ]
        
        # Number of frames
        self.nframes = len(self.frames)

        # Lenghts of frames
        self.lengths = [ j-i for i,j in zip(self.frames,self.frames[1:]+[self.nframes]) ]             

        # Make a selection
        if selection:
            self.frames  = [ j for i in selection for j in self.frames[i]  ]
            self.lengths = [ j for i in selection for j in self.lengths[i] ]

        if reverse:
            # Reverse the lists
            self.frames.reverse()
            self.lengths.reverse()

        # Get the step and time for each frame 
        if steps or times:
            self.steps, self.times = zip(*[self.file.seek(i+8) or unpack(">lf",self.file.read(8)) for i in self.frames])
            # Convert from ps to ns
            self.times = [ 0.001*i for i in self.times ]                
        else:
            self.steps = range(len(self.frames))
            self.times = range(len(self.frames))
            

# Read every file into an XTC (trajectory) object
xtcList = [ XTC(i,last=True) for i in FNM ]


# To stitch trajectories together, we need to sort the files according to the start time
if options["-stitch"]:
    xtcList.sort(key=lambda xtc: xtc.start)
    if options["-reverse"]:
        xtcList.reverse()
    for i in xtcList: print repr(i)


if options["-info"]:
    for xtc in xtcList:
        print repr(xtc)


if options["-last"]:
    with open(options["-o"].value,"wb") as out:
        for xtc in xtcList:
            out.write(xtc.frame)


if selection:
    selection = [ r(i) for j in selection for i in j.split(",") ]


if (options["-reverse"] or 
    options["-split"]   or
    options["-steps"]   or
    options["-times"]   or
    options["-stitch"]  or
    selection):

    if options["-split"]:
        # This is a pretty unsafe way of hacking the filename...
        if not "%" in options["-o"].value:
            options["-o"].value = options["-o"].value.replace(".xtc","%d.xtc")
    else:
        out  = open(options["-o"].value, "wb")

    num      = 0
    lasttime = options["-reverse"] and 1e32 or -1e32

    for xtc in xtcList:

        # Extract data for all frames
        xtc.allFrames(selection=selection, 
                      steps=options["-steps"], 
                      times=options["-times"] or options["-stitch"],
                      reverse=options["-reverse"])

        if options["-steps"]:
            print xtc.name, " ".join([str(i) for i in xtc.steps])
            continue

        if options["-times"]:
            print xtc.name, " ".join([str(i) for i in xtc.times])
            continue

        for tm, f in xtc:
            if options["-stitch"]:
                if options["-reverse"]:
                    if tm >= lasttime:
                        continue
                elif tm <= lasttime:
                    continue
            if options["-split"]:
                out = open(options["-o"].value%num,"wb")
                num += 1
            out.write(f)
            if options["-split"]:
                out.close() 
            lasttime = tm

    if not options["-split"]:
        out.close()



