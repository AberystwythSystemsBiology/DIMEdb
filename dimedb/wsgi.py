'''
  ___ ___ __  __ ___    _ _
 |   \_ _|  \/  | __|__| | |__
 | |) | || |\/| | _|/ _` | '_ \
 |___/___|_|  |_|___\__,_|_.__/

 Welcome to DIMEdb.

 Startup command: gunicorn --bind 0.0.0.0:80 wsgi:app

 Homepage URL: http://dimedb.ibers.aber.ac.uk
 Documentation URL: http://www.github.com/KeironO/DIMEdb

'''


from app import app
import logging
from logging.handlers import RotatingFileHandler

if __name__ == "__main__":
    print __doc__
    handler = RotatingFileHandler('dimedb_log.log', backupCount=1)
    handler.setLevel(logging.INFO)
    log = logging.getLogger('werkzeug')
    log.setLevel(logging.DEBUG)
    log.addHandler(handler)
    app.logger.addHandler(handler)
    app.run()