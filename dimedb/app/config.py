import os, csv

class BaseConfig(object):


    def load_settings():
        settings = {}
        fp = os.path.join(os.path.expanduser("~"), "dimedb.csv")
        with open(fp, "rb") as settings_file:
            csvreader = csv.reader(settings_file, delimiter=",")
            for row in csvreader:
                settings[row[0]] = row[1]
        return settings

    settings = load_settings()

    DEBUG = True
    BASE_DIR = os.path.abspath(os.path.dirname(__file__))
    SECRET_KEY = settings["SECRET_KEY"]

    MAINTENANCE = False

    SQLALCHEMY_DATABASE_URI = "postgresql+psycopg2://%(username)s:%(password)s@localhost/dimedb" % {
        "username" : settings["PSQL_USERNAME"],
        "password" : settings["PSQL_PASSWORD"]
    }

    SQLALCHEMY_TRACK_MODIFICATIONS = False
    LDAP_LOGIN_VIEW = 'auth.login'

    CSRF_ENABLED = True
    CSRF_SESSION_KEY = settings["CSRF_SESSION_KEY"]

    MAIL_SERVER = "smtp.googlemail.com"
    MAIL_PORT = 465
    MAIL_USE_TLS = False
    MAIL_USE_SSL = True

    SECURITY_PASSWORD_SALT = settings["SECURITY_PASSWORD_KEY"]

    MAIL_USERNAME = settings["MAIL_USERNAME"]
    MAIL_PASSWORD = settings["MAIL_PASSWORD"]

    MAIL_DEFAULT_SENDER = "dimedbemail@gmail.com"
