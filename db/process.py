import os, urllib, json

def get_periodic_table(fp="./dl-files/periodic_table.json"):
    if os.path.isfile(fp) == False:
        print "Cannot find the periodic file in %s, downloading it now." % fp
        pt = urllib.URLopener()
        pt.retrieve("https://gist.githubusercontent.com/KeironO/a2ce6d7fb7e7e10f616a51f511cb27b4/raw/71b0ad0abf92ef101a4833dd9ee0837609939f9d/gistfile1.txt", fp)
        print "File downloaded..."
        get_periodic_table()
    else:
        with open(fp, "r") as pt_fp:
            pt = json.load(pt_fp)
        return pt

if __name__ == "__main__":
    pt = get_periodic_table()
