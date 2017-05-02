from eve import Eve

from flask import render_template, abort, request, jsonify, url_for

app = Eve(__name__)
app.config.update(
    DEBUG = True
)

# Annoying webpage stuff.

@app.route("/")
def homepage():
    return render_template("main.html")

@app.route("/search/")
def search():
    return render_template("search.html")


@app.route("/help/")
def help():
    return render_template("help.html", url = request.url_root)

@app.route("/view/<string:_id>/")
def view(_id):
    try:
        return render_template("view.html" , id=_id)
    except Exception, err:
        abort(403)

@app.errorhandler(404)
def page_not_found(e):
    return render_template("./misc/404.html"), 404

# Temporary until I find a better solution.
@app.route("/gen_structure/<string:id>/")
def smiles_to_2d(id):
    from rdkit import Chem
    from rdkit.Chem import Draw
    import StringIO
    try:
        smiles = d.MetaboliteFull.objects.filter(_id=id)[0].smiles
        smiles_image = StringIO.StringIO()
        mol = Chem.MolFromSmiles(smiles)
        Draw.MolToFile(mol, fileName=smiles_image, imageType="png", size=(300, 300))
        return smiles_image.getvalue().encode("base64")
    except Exception, err:
        abort(404)
