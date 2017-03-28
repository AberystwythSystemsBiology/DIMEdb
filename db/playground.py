from bioservices import *
import urllib2

compound = "C00089"

url = "http://www.genome.jp/dbget-bin/get_linkdb?-t+genes+cpd:"+compound

response = urllib2.urlopen(url)
html = response.read()

print html