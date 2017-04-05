function getBaseURL () {
   return location.protocol + "//" + location.hostname +
      (location.port && ":" + location.port) + "/";
}

// Generalised Results Generator
function render_search_results(table_id, api_url, length) {
    $('#'+table_id+'').delay(100).DataTable({
                    "destroy": true,
                    "ajax": getBaseURL()+api_url,
                    "columns": [
                        {
                            "title": "Metabolite Name",
                            "width": "60%",
                            "data": "name",
                            "render": function (data, type, row) {
                                return data

                            }
                        },
                        {
                            "title": "Molecular Formula",
                            "width": "10%",
                            "className": "dt-center",
                            "data": "molecular_formula",
                            "render": function (data, type, row) {
                                return '<p style="text-align:center">'+data.replace(/([0-9]+)/g, '<sub>$1</sub>')+'</p>';
                            }
                        },
                        {
                            "title": "Neutral Mass (m/z)",
                            "data": "accurate_mass",
                            "className": "dt-center",
                            "width": "10%",
                            "render": function (data, type, row) {
                                return data.toFixed(4);
                            }
                        },
                        {
                            "title": "Actions",
                            "data": "id",
                            "className": "dt-center",
                            "width": "15%",
                            "render": function (data, type, row) {
                                var view_url = getBaseURL() + "view/" + data;
                                return "<div id='tester' class='btn-toolbar text-center'>" +
                                    "<a href='" + view_url + "' target='_blank'>" +
                                    "<button class='btn btn-sm btn-primary' id='view_button'>View</button>" +
                                    "</a>" +
                                    "<button id='view_card' class='btn btn-sm btn-danger' name='"+data+"'><i class='glyphicon glyphicon-book'></i></button>" +
                                    "<button id='add_to_clipboard' class='btn btn-sm btn-success' name='"+data+"'><i class='glyphicon glyphicon-plus'></i></button></div>"
                            }
                        }
                    ],
                    "searching": false,
                    //"bSort" : false,
                    "lengthChange": false,
                    "pageLength": length
                });
}

function clear_card() {
    $("#card_metabolite_name").empty();
    $("#molecular_formula").empty();
    $("#neutral_mass").empty();
    $("#positive_adduct").empty();
    $("#negative_adduct").empty();
    $("#origins").empty();
    $("#biofluids").empty();
    $("#tissues").empty();
}

function generate_card(id) {
    var current_url = getBaseURL();
    $.getJSON(current_url + "api/metabolite/?id=" +id, function(json){
        var metabolite = json["data"][0];
        clear_card();

        $("#card_metabolite_name").html(metabolite["name"]);
        $("#molecular_formula").html(metabolite["molecular_formula"].replace(/([0-9]+)/g, '<sub>$1</sub>'));

        for (result in metabolite["adducts"]) {
            if (result == "neutral") {
                $("#neutral_mass").html(metabolite["adducts"]["neutral"]["peaks"][0]["accurate_mass"]);
            }

            else {
                for (adduct_idx in metabolite["adducts"][result]["peaks"]) {
                    var adduct = metabolite["adducts"][result]["peaks"][adduct_idx];
                    $("#" + result + "_adduct").append("<li class='list-group-item'><b>" + adduct["type"] + ":</b> " + adduct["accurate_mass"].toFixed(4)+"</li>");
                }
            }
        }

        if (metabolite["origins"] == null) {
            $("#origins").append("<i class='text-primary'>None</i>");
        }
        else {
            for (indx in metabolite["origins"]) {
                $("#origins").append(metabolite["origins"][indx] + "; ");
            }
        }

        if (metabolite["biofluid_locations"] == null) {
            $("#biofluids").append("<i class='text-primary'>None</i>");
        }
        else {
            for (indx in metabolite["biofluid_locations"]) {
                $("#biofluids").append(metabolite["biofluid_locations"][indx] + "; ");
            }
        }

        if (metabolite["tissue_locations"] == null) {
            $("#tissues").append("<i class='text-primary'>None</i>");
        }
        else {
            for (indx in metabolite["tissue_locations"]) {
                $("#tissues").append(metabolite["tissue_locations"][indx] + "; ");
            }
        }

        $("#viewfromcard").attr("href", current_url+"view/"+id);

    });
    $("#metabolite_card").modal("toggle");
}