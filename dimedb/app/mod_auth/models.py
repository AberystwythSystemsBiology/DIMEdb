from app import db, bcrypt
from sqlalchemy.ext.hybrid import hybrid_property
from hashlib import md5

class Base(db.Model):
    __abstract__ = True

    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    creation_date = db.Column(db.DateTime, default=db.func.current_timestamp())
    modification_date = db.Column(db.DateTime, default=db.func.current_timestamp(), onupdate=db.func.current_timestamp())

class UserType(db.Enum):
    ADMIN = "Administrator"
    RESEARCHER = "Researcher"
    OTHER = "Other"


class User(Base):
    __tablename__ = "auth_user"
    user_name = db.Column(db.String(64), unique=True, nullable=False)
    _password = db.Column(db.String(128), nullable=False)
    email_address = db.Column(db.String(128), unique=True)

    affiliation = db.Column(db.String(128))
    country = db.Column(db.String(128))

    first_name = db.Column(db.String(128), nullable=False)
    last_name = db.Column(db.String(128), nullable=False)
    user_type = db.Column(db.String, nullable=False)

    @hybrid_property
    def password(self):
        return self._password

    @password.setter
    def _set_password(self, plaintext):
        self._password = bcrypt.generate_password_hash(plaintext)

    def is_correct_password(self, plaintext):
        if bcrypt.check_password_hash(self._password, plaintext):
            return True
        else:
            return False

    def is_authenticated(self):
        return True

    def is_active(self):
        return True

    def is_anonymous(self):
        return False

    def get_id(self):
        return unicode(self.id)

    def avatar(self, size):
        return 'https://www.gravatar.com/avatar/%s?d=mm&s=%d' % (
        md5(self.email_address.encode('utf-8')).hexdigest(), size)

    def __repr__(self):
        return self