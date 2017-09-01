from flask_wtf import FlaskForm

from wtforms import StringField, RadioField, BooleanField
from wtforms.validators import DataRequired, EqualTo, Length

class CreateTableForm(FlaskForm):
    title = StringField("Title", [DataRequired("You need to provide a title."), Length(6)])
    description = StringField("Description", [Length(max=6000)])
    species = StringField("Species", [Length(max=512)])

class EditTableForm(FlaskForm):
    title = StringField("Title", [DataRequired("You need to provide a title."), Length(6)])
    description = StringField("Description", [Length(max=6000)])
    species = StringField("Species", [Length(max=512)])

class NewPublicationForm(FlaskForm):
    doi = StringField("DOI")
    pubmed_id = StringField("Pubmed ID")
    publication_title = StringField("Publication Title", [DataRequired("You need to provide a title.")])
    author_list = StringField("Author List", [DataRequired("You need to provide a author list.")])
