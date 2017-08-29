from flask.ext.script import Manager
from flask_migrate import Migrate, MigrateCommand

from app import app, db
from app.mod_auth.models import User

app.config_from_object(app.config)

migrate = Migrate(app, db)
manager = Manager(app)

manager.add_command('db', MigrateCommand)

@manager.command
def create_db():
    db.create_all()

@manager.command
def drop_db():
    db.drop_all()

if __name__ == "__main__":
    manager.run()