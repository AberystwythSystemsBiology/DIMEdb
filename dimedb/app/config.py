import os

class BaseConfig(object):

    config_dict = {
        "psql_un" : os.environ["DIMEDB_PSQL_USERNAME"],
        "psql_p" : os.environ["DIMEDB_PSQL_PASSWORD"]
    }

    DEBUG = True
    BASE_DIR = os.path.abspath(os.path.dirname(__file__))
    SECRET_KEY = os.environ["SECRET_KEY"]

    MAINTENANCE = False

    SQLALCHEMY_DATABASE_URI = "postgresql+psycopg2://%(psql_un)s:%(psql_p)s@localhost/dimedb" % config_dict
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    LDAP_LOGIN_VIEW = 'auth.login'

    CSRF_ENABLED = True
    CSRF_SESSION_KEY = os.environ["CSRF_SESSION_KEY"]

    MAIL_SERVER = "smtp.googlemail.com"
    MAIL_PORT = 465
    MAIL_USE_TLS = False
    MAIL_USE_SSL = True

    SECURITY_PASSWORD_SALT = os.environ["SECURITY_PASSWORD_KEY"]

    MAIL_USERNAME = os.environ["DIMEDB_MAIL_USERNAME"]
    MAIL_PASSWORD = os.environ["DIMEDB_MAIL_PASSWORD"]

    MAIL_DEFAULT_SENDER = "dimedbemail@gmail.com"