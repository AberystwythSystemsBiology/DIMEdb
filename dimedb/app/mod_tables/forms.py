from flask_wtf import FlaskForm

from wtforms import StringField, RadioField, BooleanField
from wtforms.validators import DataRequired, EqualTo, Length

class CreateTableForm(FlaskForm):
    title = StringField("Title", [DataRequired("You need to provide a title."), Length(6)])
    description = StringField("Description", [Length(max=6000)])
    species = StringField("Species", [Length(max=512)])
    doi = StringField("DOI")
    public = BooleanField("Public")

class EditTableForm(FlaskForm):
    title = StringField("Title", [DataRequired("You need to provide a title."), Length(6)])
    description = StringField("Description", [Length(max=6000)])
    species = StringField("Species", [Length(max=512)])
    doi = StringField("DOI")