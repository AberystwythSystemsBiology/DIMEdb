from app import db

class Base(db.Model):
    __abstract__ = True

    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    creation_date = db.Column(db.DateTime, default=db.func.current_timestamp())
    modification_date = db.Column(db.DateTime, default=db.func.current_timestamp(), onupdate=db.func.current_timestamp())

class MetaboliteTables(Base):
    __tablename__ = "metabolite_tables"
    owner_id = db.Column(db.Integer, db.ForeignKey("auth_user.id"), nullable=False)
    title = db.Column(db.String, nullable=False)
    description = db.Column(db.String, nullable=True)
    doi = db.Column(db.String, nullable=True)
    public = db.Column(db.Boolean, default=False, nullable=False)

class Metabolite(Base):
    __tablename__ = "metabolite"
    table_id = db.Column(db.Integer, db.ForeignKey("metabolite_tables.id"), nullable=False)
    inchikey = db.Column(db.String, nullable=False)
    comment = db.Column(db.String, nullable=True)