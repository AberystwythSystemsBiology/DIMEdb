import os

class BaseConfig(object):

    config_dict = {
        "psql_un" : os.environ["DIMEDB_PSQL_USERNAME"],
        "psql_p" : os.environ["DIMEDB_PSQL_PASSWORD"],
        "secret_k" : os.environ["SECRET_KEY"],
        "csrf_s_k" : os.environ["CSRF_SESSION_KEY"],
        "security_p_s" : os.environ["SECURITY_PASSWORD_KEY"]
    }

    DEBUG = True
    BASE_DIR = os.path.abspath(os.path.dirname(__file__))
    SECRET_KEY = config_dict["secret_k"]

    MAINTENANCE = False

    SQLALCHEMY_DATABASE_URI = "postgresql+psycopg2://%(psql_un)s:%(psql_p)s@localhost/dimedb" % config_dict
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    LDAP_LOGIN_VIEW = 'auth.login'

    CSRF_ENABLED = True
    CSRF_SESSION_KEY = config_dict["csrf_s_k"]

    MAIL_SERVER = "smtp.googlemail.com"
    MAIL_PORT = 465
    MAIL_USE_TLS = False
    MAIL_USE_SSL = True

    SECURITY_PASSWORD_SALT = config_dict["security_p_s"]

    MAIL_USERNAME = os.environ["DIMEDB_MAIL_USERNAME"]
    MAIL_PASSWORD = os.environ["DIMEDB_MAIL_PASSWORD"]

    MAIL_DEFAULT_SENDER = "dimedbemail@gmail.com"