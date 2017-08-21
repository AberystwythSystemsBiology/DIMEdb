from functools import wraps
from flask import flash, redirect, url_for, request, abort
from flask_login import current_user

def check_confirmed(func):
    @wraps(func)
    def decorated_function(*args, **kwargs):
        if current_user.is_confirmed() != True:
            flash("Awaiting account confirmation")
            return request.referrer
        return func(*args, **kwargs)
    return decorated_function