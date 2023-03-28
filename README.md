# metabolitephotochem
Welcome! This is the repository for the code that I used to assess metabolite degradation rates under a solar simulator. 

## A Brief Set of Notes
Hello there! If you're viewing this repository, well, it's a work in progress. I'm still in the process of making it self-consistent, as currently the code here reflects the state of analysis while I'm still toying with overly ambitious ideas. 
- There may be data here that I use for things that won't be in the final paper, as well as accompanying scripts and functions.
- Since this is a code repository, the larger data files are _not here._ I'm keeping them locally in a very similar file structure, and should you need the whole thing, ask!
- My "versioning" for actual `.mat`-based data files is often date-based, meaning the references to specific files and directories in the code will need to be checked if you run things. I've included an example dataset that will work from the `processing.m` script onwards.

## This Repository Has a Dependency...
...but I've included it. If the idea of a bunch of color palettes based on album covers appeals to you, you can download it separately and save it to your MATLAB home directory [here.](https://github.com/germo006/NoahMaps)

## Down to Business

### If you're curious enough to work through this from the raw (or raw-ish) data...
...start from the beginning. 
1. [Contact me](ngermolus@whoi.edu) and I can get you data I processed in Skyline, or you can download the .RAW files yourself from the accompanying MetaboLights repository [here](www.ebi.ac.uk/metabolights/MTBLS7513) and pick peaks! Go nuts.
2. If you've started from there, you will end up with a list of light and heavy peak areas. I've been using Skyline a while, but expect a couple days of monotonous work here. Unfortunately, you'll need to use [Skyline](https://skyline.ms/project/home/begin.view) nomenclature for the columns, _and_ have a sequence file (basically the list of samples ran on the Orbitrap with some metadata). That can be found in the [datasets folder](../blob/master/datasets).
3. You'll _also_ need the two "transition lists" containing information about each molecule, which is in that [datasets folder](../blob/master/datasets). There's one for positive and negative ion modes.
4. Finally, you should be able to run the routine `considerSkyline.m` from its wrapper, `readSkyline.m`. You'll have to check filepaths to make sure they work with your local machine, but this will spit out actual concentration data in pM units.

### If you're using my example dataset or did the above...
...Check out `processing.m`. This script is long, but it serves as the main wrapper. Turn on all the `for` loops, and you should be able to generate every single figure I included in the paper, as well as many that I did not.
- There is a whole section or two (and a several functions called here) that are not in the paper. They concern quantum yields. That was cut from the paper for a reason, and you should disregard it. It is not useful.