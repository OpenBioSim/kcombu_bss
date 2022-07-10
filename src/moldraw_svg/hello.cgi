#!/usr/bin/env python

import sys
import cgi
import os
from datetime import datetime

now = datetime.now()
date_string = now.strftime("%Y/%m/%d %H:%M:%S")


print "Content-type: text/html\n";
print "<HTML><BODY><PRE>";
print "HELLO !! This is python cgi-bin(%s)."%(date_string);
form = cgi.FieldStorage()
optdic = {}
for key in (form.keys()):
  optdic[key] = form[key].value
  print "key:'%s' value:'%s'"%(key,form[key].value)


print "</PRE></BODY></HTML>"
