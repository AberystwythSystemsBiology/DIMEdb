// Need to write a API handler to retrieve all data.

function render_view(base_url, id) {
    var url = base_url+"api/metabolites/?id__exact="+id;
    $.getJSON(url, function (data) {
        var metabolite = data["data"][0];
        // Change the document title.
        $(document).prop('title', metabolite["name"] + " : Metabolite DataBase");
        $("#name").html(metabolite["name"]);
        $("#mf").html(metabolite["molecular_formula"].replace(/([0-9]+)/g, '<sub>$1</sub>'));
        $("#am").html(metabolite["accurate_mass"]);
        $("#o").html(metabolite["origins"]);
        $("#s").html(metabolite["smiles"]);
    });



}