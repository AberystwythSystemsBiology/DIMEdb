from flask_wtf import FlaskForm

from wtforms import StringField, RadioField, TextField
from wtforms.validators import DataRequired, EqualTo, Length

class CreateTableForm(FlaskForm):
    title = StringField("Title", [DataRequired("Forgot your user name?"), Length(6)])
    description = StringField("Description", [Length(max=6000)])
    species = StringField("Species", [Length(max=512)])
    doi = StringField("DOI")
    public = RadioField("Public", choices=[("True", True), ("False", False)])