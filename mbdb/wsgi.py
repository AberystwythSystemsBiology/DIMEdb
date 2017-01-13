#gunicorn --bind 0.0.0.0:80 wsgi:app

from app import app

if __name__ == "__main__":
    app.run()
