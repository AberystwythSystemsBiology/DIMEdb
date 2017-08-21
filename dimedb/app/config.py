import os

class BaseConfig(object):
    DEBUG = True
    BASE_DIR = os.path.abspath(os.path.dirname(__file__))
    SECRET_KEY = "secret"

    SQLALCHEMY_DATABASE_URI = "postgresql+psycopg2://keo7:password@localhost/dimedb"
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    LDAP_LOGIN_VIEW = 'auth.login'
    CSRF_ENABLED = True
    CSRF_SESSION_KEY = "secret"

    MAIL_SERVER = "smtp.googlemail.com"
    MAIL_PORT = 465
    MAIL_USE_TLS = False
    MAIL_USE_SSL = True

    SECURITY_PASSWORD_SALT = "salty"

    MAIL_USERNAME = os.environ["DIMEDB_MAIL_USERNAME"]
    MAIL_PASSWORD = os.environ["DIMEDB_MAIL_PASSWORD"]

    MAIL_DEFAULT_SENDER = "dimedbemail@gmail.com"