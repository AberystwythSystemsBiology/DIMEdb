from itsdangerous import URLSafeSerializer

from app import app

def generate_confirmation_token(email):
    serialiser = URLSafeSerializer(app.config["SECRET_KEY"])
    return serialiser.dumps(email, salt=app.config["SECURITY_PASSWORD_SALT"])

def confirm_token(token, expiration=3600):
    serialiser = URLSafeSerializer(app.config["SECRET_KEY"])
    try:
        email = serialiser.loads(token, salt=app.config["SECURITY_PASSWORD_SALT"], max_age=expiration)
    except:
        return False
    return email