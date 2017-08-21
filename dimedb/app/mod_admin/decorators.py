from functools import wraps
from flask import flash, redirect, url_for, request, abort
from flask_login import current_user

def check_admin(func):
    @wraps(func)
    def decorated_function(*args, **kwargs):
        if current_user.is_admin() != True:
            abort(403)
        return func(*args, **kwargs)
    return decorated_function