'''
  ___ ___ __  __ ___    _ _
 |   \_ _|  \/  | __|__| | |__
 | |) | || |\/| | _|/ _` | '_ \
 |___/___|_|  |_|___\__,_|_.__/

 Welcome to DIMEdb.


 Homepage URL: http://dimedb.ibers.aber.ac.uk
 Documentation URL: http://www.github.com/KeironO/DIMEdb

'''


from dimedb import app
import logging
from logging.handlers import RotatingFileHandler

if __name__ == "__main__":
    print __doc__
    app.run()
