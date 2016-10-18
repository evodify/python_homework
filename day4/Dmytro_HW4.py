#!/usr/bin/python2

import os, warnings, webbrowser

# Home work assignment 4

class Song(object):
  """It's a class that describes a person"""
  def __init__(self,title,artist,duration): # attributes
    self.title = title
    self.artist = artist
    self.duration = duration
    try:
      self.duration = int(duration)
    except:
      #warnings.warn('Durations of songs that have non-numerical duration are set to 0.')
      print('Warning: "%s" song has non-numerical duration of "%s". It is set to 0.' % (self.title, self.duration))
      self.duration = 0
    if self.duration < 0:
      #raise Exception('"%s" song has negative duration of "%i"' % (self.title, self.duration))
      print('Warning: "%s" song has negative duration of "%i". It is set to 0.' % (self.title, self.duration))
      self.duration = 0
  def pretty_duration(self): # method
    """Returns a nice string describing the duration. For instance if the duration is 3611, this methods takes no input and returns "01 hours 00 minutes 11 seconds" as a string."""
    hour = self.duration/360
    minut = (self.duration-hour*360)/60
    sec = self.duration-hour*360-minut*60
    return "%i hours %i minutes %i seconds" % (hour, minut, sec)
  def play(self):
    """Automatically opens a webpage on your computer with a youtube search for the title."""
    artistTitle = self.artist+' '+self.title
    quaryArtistTitle = "+".join(artistTitle.split())
    url = "https://www.youtube.com/results?search_query="+quaryArtistTitle # create an url
    webbrowser.open(url, new = 2) # open in a browser in a new tab

in_path = "lulu_mix_16.csv"
# Check that the path is valid #
if not os.path.exists(in_path):
  raise Exception("No file at '%s'." % in_path)

with open(in_path) as datafile:
  header_line = datafile.readline()
  header_words = header_line.split(',')
  songs = []
  for line in datafile:
    words = line.split(',')
    songs.append(Song(title = words[0],artist = words[1],duration = words[2]))

for s in songs: print s.duration
for s in songs: print s.artist
for s in songs: print s.pretty_duration()
print sum(s.duration for s in songs), "seconds in total"
songs[6].play()
