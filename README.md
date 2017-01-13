# bisc484

## Synopsis

Companion web app for UD's BISC 484: Computer Based Genetics Lab. Code implements a flask based web app.

## Installation

bioApp/

Entire directory is placed in /srv. User:Group should be changed to uwsgi:nginx for entire directory. For SELinux, context should be set to system\_u:object\_r:unlabeled\_t:s0 . As a last resort, httpd context can be set to permissive: semanage permissive -a httpd\_t .


