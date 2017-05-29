#!/usr/bin/env python
"""
This script crawls several news sources to construct an encyclopedia of RNA objects. 
"""

import feedparser
from urllib2 import urlopen
from pyrna.websites import Pubmed, RNA, RNABiology, NAR

def find_non_business_words():
    """
    This function parses feeds of several websites to extract the english words that are not related to any business domain.
    """
    all_words = None
    for url in ["http://ieeexplore.ieee.org/rss/TOC7.XML", "http://ieeexplore.ieee.org/rss/TOC85.XML", "http://ieeexplore.ieee.org/rss/TOC8.XML", "http://ieeexplore.ieee.org/rss/TOC6221020.XML"]:
        feed = feedparser.parse(url)
        for entry in feed['entries']:
            new_words = set([word.lower() for word in entry.description.split()])
            if not all_words:
                all_words =  new_words
            else:
                all_words = all_words.union(new_words)            
    return all_words

def crawl():
    #l = list(non_domain_specific_words)
    #l.sort()
    #print l
    return Pubmed().parse().union(RNA().parse()).union(RNABiology().parse()).union(NAR().parse())

if __name__ == '__main__':
    non_business_words = find_non_business_words()
    business_words = crawl().difference(non_business_words)
    print business_words
    print len(business_words)