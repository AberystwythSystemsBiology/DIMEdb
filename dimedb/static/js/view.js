function getBaseURL () {
   return location.protocol + "//" + location.hostname +
      (location.port && ":" + location.port) + "/";
}

function get_metabolite(id) {
    var api_url = encodeURI(getBaseURL()+'api/metabolites/?where={"_id" : "'+id+'"}');

    var json = (function () {
        var json = null;
        $.ajax({
            'async': false,
            'global': false,
            'url': api_url,
            'dataType': "json",
            'success': function (data) {
                json = data;
            }
        });
        return json;
    })();
    return json["_items"][0]
}

function basic_information(metabolite) {

    if (metabolite["synonyms"].length != 0) {
        var synonyms = "";
        for (i = 0; i < metabolite["synonyms"].length; i++) {
            synonyms += metabolite["synonyms"][i] + "; ";
        }
    }
    else {
        synonyms = "Not available";

    }

    $("#synonyms").html(synonyms);

     $("#met_name").html(metabolite["name"]);
     $("#met_name_modal").html(metabolite["name"]);
     $("#molecular_formula").html(metabolite["molecular_formula"].replace(/([0-9]+)/g, '<sub>$1</sub>'));
     $("#accurate_mass").html(metabolite["accurate_mass"].toFixed(6));
     $("#neutral_mass").html(metabolite["adducts"]["neutral"][0]["accurate_mass"].toFixed(6));

     console.log(metabolite["name"]);
}

function render_view(id) {
    var metabolite = get_metabolite(id);
    basic_information(metabolite)
}