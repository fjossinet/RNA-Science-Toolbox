import feedparser

class Twitter:
	def __init__(self, profile_name):
		self.profile_name = profile_name

	def parse(self):
		pass

class Website:
	
	def __init__(self, web_pages = {}, feeds = {}):
		self.web_pages = web_pages
		self.feeds = feeds

	def parse(self):
		all_words = None
		#rss channels parsed with feedparser
		for name,url in self.feeds.items():
			feed = feedparser.parse(url)
			for entry in feed['entries']:
				new_words = set([word.lower() for word in entry.title.split()])
				if not all_words:
					all_words =  new_words
				else:
					all_words = all_words.union(new_words)
		return all_words

class NAR(Website):

	def __init__(self):
		Website.__init__(self, web_pages = {}, feeds = {"recent issues":"http://nar.oxfordjournals.org/rss/recent.xml", "advance access": "http://nar.oxfordjournals.org/rss/ahead.xml"})

	def parse(self):
		return Website.parse(self)

class Pubmed(Website):

	def __init__(self):
		Website.__init__(self, web_pages = {}, feeds = {"ncRNA":"http://www.ncbi.nlm.nih.gov/entrez/eutils/erss.cgi?rss_guid=1BwDT-LAhNxOxdpRxspUhu7jKysyoW_xbFNHPji6mc7Et8lcgM"})

	def parse(self):
		return Website.parse(self)

class RNA(Website):
	def __init__(self):
		Website.__init__(self, web_pages = {}, feeds = {"current issue": "http://rnajournal.cshlp.org/rss/current.xml", "in advance": "http://rnajournal.cshlp.org/rss/ahead.xml"})

	def parse(self):
		return Website.parse(self)

class RNABiology(Website):
	def __init__(self):
		Website.__init__(self, web_pages = {}, feeds = {"current issue": "http://www.landesbioscience.com/rss/journals/rnabiology/current/"})

	def parse(self):
		return Website.parse(self)





