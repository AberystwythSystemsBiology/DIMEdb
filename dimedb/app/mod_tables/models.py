from dimedb import db

class Base(db.Model):
    __abstract__ = True

    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    creation_date = db.Column(db.DateTime, default=db.func.current_timestamp())
    modification_date = db.Column(db.DateTime, default=db.func.current_timestamp(), onupdate=db.func.current_timestamp())

class MetaboliteTable(Base):
    __tablename__ = "metabolite_table"
    owner_id = db.Column(db.Integer, db.ForeignKey("auth_user.id"), nullable=False)
    title = db.Column(db.String(2048), nullable=False)
    description = db.Column(db.String(2048), nullable=True)
    published = db.Column(db.Boolean, default=False, nullable=False)
    removed = db.Column(db.Boolean, default=False, nullable=False)
    species = db.Column(db.String, nullable=True)
    locked = db.Column(db.Boolean, nullable=False, default=False)

class MetaboliteTablePublication(Base):
    __tablename__ = "metabolite_table_publication"
    table_id = db.Column(db.Integer, db.ForeignKey("metabolite_table.id"), nullable=False)
    pubmed_id = db.Column(db.String(30))
    doi = db.Column(db.String(100))
    author_list = db.Column(db.String(100))
    publication_title = db.Column(db.String(2048))
    publication_status = db.Column(db.String())

class Metabolite(Base):
    __tablename__ = "metabolite"
    table_id = db.Column(db.Integer, db.ForeignKey("metabolite_table.id"), nullable=False)
    inchikey = db.Column(db.String, nullable=False)
    polarity = db.Column(db.String, nullable=True)
    adduct = db.Column(db.String, nullable=True)
    comments = db.Column(db.String(2048))
