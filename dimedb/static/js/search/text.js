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
    var api_url = getBaseURL()+'api/metabolites/?where={ "$and" : [';
    /*
    if (values["Name"] != null) {
        console.log(values["Name"]);
        api_url += '"$or" : [';
        api_url += '{"Identification Information.Name" :';
        api_url += generate_regex(values["Name"]) + "},";
        api_url += '{"Identification Information.IUPAC Name" :';
        api_url += generate_regex(values["Name"]) + "},";
        api_url += '{"Identification Information.Systematic Name" :';
        api_url += generate_regex(values["Name"]) + "},";
        api_url += '{"Identification Information.Synonyms" : {"$in" : [';
        api_url += generate_regex(values["Name"]) + "]}}";
        api_url += "]"
    }
    */

    if (values["Molecular Formula"] != null) {
        api_url = check_if_comma(api_url);
        api_url += '{"Identification Information.Molecular Formula" : ';
        api_url += is_equals_to(values["Molecular Formula"]);
        api_url += '}';
    }

    if (values["Molecular Weight"] != null) {
        var tolerance = parseFloat($("#search_tolerance").val());
        var lte = parseFloat(values["Molecular Weight"]) + tolerance;
        var gte = parseFloat(values["Molecular Weight"]) - tolerance;
        api_url = check_if_comma(api_url);
        api_url += '{"Physiochemical Properties.Molecular Weight" : ';
        api_url += '{"$lte" : ' + lte + ', "$gte" : ' + gte + '}';
        api_url += '}';
    }

    if (values["External Source"] != null) {
        var source_name = $("#source_name").val();
        api_url = check_if_comma(api_url);
        api_url += '{"External Sources.'+source_name+'" : "';
        api_url += values["External Source"] + '"}';

    }

    if (values["InChI"] != null) {
        api_url = check_if_comma(api_url)
        api_url += '{"Identification Information.InChI" : ';
        api_url += is_equals_to(values["InChI"]);
        api_url += '}';
    }

    if (values["InChI Key"] != null) {
        api_url = check_if_comma(api_url)
        api_url += '{"_id" : ';
        api_url += is_equals_to(values["InChI Key"]);
        api_url += '}';
    }

    if (values["SMILES"] != null) {
        api_url = check_if_comma(api_url)
        api_url += '{"Identification Information.SMILES" : ';
        api_url += is_equals_to(values["SMILES"]);
        api_url += '}';
    }

    api_url += "]}";

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
    console.log(values);
    if (checkProperties(values) == false) {
        var api_url = generate_api_url(values);
        results_table(api_url);

    }
    $("#advanced_search_results").fadeIn(200);
});
