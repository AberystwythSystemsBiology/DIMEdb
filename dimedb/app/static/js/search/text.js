function text_to_null(text) {
    if (text == "") {
        return null
    }
    else {
        return $.trim(text)
    }
}

function get_values() {
    return {
        "Name" : text_to_null($("#name").val()),
        "Molecular Formula" : text_to_null($("#molecular_formula").val()),
        "Molecular Weight" : text_to_null($("#molecular_weight").val()),
        "External Source" : text_to_null($("#source_id").val()),
        "InChI" : text_to_null($("#inchi").val()),
        "InChI Key" : text_to_null($("#inchikey").val()),
        "SMILES" : text_to_null($("#smiles").val())
    }
}

function contains_regex(s) {
    return '{"$regex":"'+s+'.*"}';
}

function is_equals_to(s) {
    return '{"$regex":"^(?i)'+s+'.*"}';
}

function checkProperties(obj) {
    for (var key in obj) {
        if (obj[key] !== null && obj[key] != "")
            return false;
    }
    return true;
}

function check_if_comma(s) {
    if (s.slice(-1) != "}") {
        return s;
    }
    else {
        s += ",";
        return s;
    }
}

function generate_api_url(values) {
    var api_url = getBaseURL() + 'api/metabolites?where={"$and" : [';

    if (values["Name"] != null) {
        var name_string = '{"$or" : [{"Identification Information.Name" : {"$regex" : ".*%s.*"}},' +
            '{"Identification Information.Synonyms" : {"$regex" : ".*%s.*"}},' +
            '{"Identification Information.IUPAC Name" : {"$regex" : ".*%s.*"}}]}';
        name_string = name_string.replace(/%s/g, values["Name"]);
        if (api_url.slice(-1) == "}") {
            api_url += ", "
        }

        api_url += name_string
    }


    if (values["Molecular Formula"] != null) {
        var mol_form = '{"Identification Information.Molecular Formula" : {"$regex" : ".*%s.*"}}';
        mol_form = mol_form.replace(/%s/g, values["Molecular Formula"]);

        if (api_url.slice(-1) == "}") {
            api_url += ", "
        }

        api_url += mol_form

    }

    if (values["SMILES"] != null) {
        var smiles = '{"Identification Information.SMILES" : {"$regex" : ".*%s.*"}}';
        smiles = smiles.replace(/%s/g, values["SMILES"])

        if (api_url.slice(-1) == "}") {
            api_url += ", "
        }

        api_url += smiles
    }

    if (values["InChI"] != null) {
        var inchi = '{"Identification Information.InChI" : {"$regex" : ".*%s.*"}}';
        inchi = inchi.replace(/%s/g, values["InChI"])

        if (api_url.slice(-1) == "}") {
            api_url += ", "
        }

        api_url += inchi
    }

    if (values["InChI Key"] != null) {
        var inchikey = '{"_id" : {"$regex" : ".*%s.*"}}';
        inchikey = inchikey.replace(/%s/g, values["InChI Key"])

        if (api_url.slice(-1) == "}") {
            api_url += ", "
        }

        api_url += inchikey
    }

    if (values["External Source"] != null) {
        var source = '{"External Sources.%source" : "%s"}}';
        source = source.replace(/%source/g, $("#source_name").val());
        source = source.replace(/%s/g, values["External Source"]);

        if (api_url.slice(-1) == "}") {
            api_url += ", "
        }

        api_url += source
    }

    if (values["Molecular Weight"] != null) {
        var tolerance = parseFloat($("#search_tolerance").val());
        var w = parseFloat(values["Molecular Weight"]);
        console.log(tolerance);

        var weight = '{"Physicochemical Properties.Molecular Weight" : {' +
            '"$gte" : %gtev, "$lte" : %ltev}}';

        weight = weight.replace(/%gtev/g, String(w - tolerance));
        weight = weight.replace(/%ltev/g, String(w + tolerance));

        if (api_url.slice(-1) == "}") {
            api_url += ", "
        }

        api_url += weight
    }


    api_url += ']}';

    console.log(api_url);

    return api_url;

}

function results_table(api_url) {
    $('#search_results_table').DataTable({
            "destroy": true,
            "searching": false,
            "lengthChange": false,
            "pageLength": 10,
            "ajax": {
                "url": encodeURI(api_url),
                "dataSrc": "_items"
            },
            "columns": [
                {
                    "title" : "Molecular Structure",
                    "width": "10%",
                    "data" : "_id",
                    "render" : function(data,type,row) {
                        return "<img src='" + getBaseURL() + "view/structure/" + row._id + "' class='img-responsive img-circle'>"
                    }
                },
                {
                    "title": "Metabolite Name",
                    "width": "40%",
                    "data": "Identification Information.Name",
                    "render": function (data, type, row) {
                        return "<a href='" + getBaseURL() + "view/" + row._id + "' target='_blank'>" + data + "</a>"
                    }
                },
                {
                    "title": "Molecular Formula",
                    "width": "10%",
                    "className": "dt-center",
                    "data": "Identification Information.Molecular Formula",
                    "render": function (data, type, row) {
                        return data.replace(/([0-9]+)/g, '<sub>$1</sub>');
                    }
                },
                {
                    "title": "Molecular Weight (g/mol)",
                    "data": "Physicochemical Properties",
                    "className": "dt-center",
                    "width": "10%",
                    "render": function (data, type, row) {
                        return data["Molecular Weight"].toFixed(3);
                    }
                }
            ],
        });
}


$("#submit_search").click(function () {
    var values = get_values();
    if (checkProperties(values) == false) {
        $("#advanced_search_results").fadeOut(200);
        var api_url = generate_api_url(values);
        results_table(api_url);

        $("#advanced_search_results").fadeIn(200);
    }

});

$("#show_filter").click(function () {
});
