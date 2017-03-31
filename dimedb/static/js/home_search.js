$(document).ready(function () {
    $(".navbar img").hide();
    $('#home-page-search-results').hide();
    $("#search_button").click(function (e) {
        var query = $("#search_text").val();
        if (query != "") {
            render_search_results(query);
            $("#home-page-welcome-search").fadeOut();
            $("#search_text_results").val(query);
            $("#home-page-search-results").delay(460).fadeIn();
        }
        else {
            alert("You don't seem to have entered anything?");
        }
    });

    $('#search_text').keypress(function (e) {
        if (e.which == 13) {
            $('#search_button').click();
        }
    });

    $("#search_text_button").click(function (e) {
        var query = $("#search_text_results").val();
        if (query != "") {
            $("#search_results").fadeOut();
            render_search_results(query);
            $("#search_results").delay(460).fadeIn();

        }
        else {
            alert("You don't seem to have entered anything?");
        }
    });

    $('#search_text_results').keypress(function (e) {
        if (e.which == 13) {
            $('#search_text_button').click();
        }
    });

    $('#home-page-search-results').on('click', "#view_card", function() {
        $("#card_table td").html("");
        generate_card($(this).attr("name"))
    });

    $('#home-page-search-results').on('click', "#add_to_clipboard", function() {
        clipboard_hook($(this).attr("name"));
    });

});

function generate_card(id) {
    var current_url = window.location.href;
    $.getJSON(current_url + "api/metabolite/?id=" +id, function(json){
        var metabolite = json["data"][0];
        console.log(metabolite);
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
    // Get the JSON data
    $("#card").modal("toggle");
}


function render_search_results(query) {
    var current_url = window.location.href;
    var url = current_url + "api/metabolites/?name__icontains=" + String(query);
    $('#search_results').delay(100).DataTable({
                    "destroy": true,
                    "ajax": url,
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
                            "title": "",
                            "data": "id",
                            "className": "dt-center",
                            "width": "15%",
                            "render": function (data, type, row) {
                                var view_url = current_url + "view/" + data;
                                return "<div id='tester' class='btn-toolbar'>" +
                                    "<a href='" + view_url + "' target='_blank'>" +
                                    "<button class='btn btn-sm btn-primary' id='view_button'>View</button>" +
                                    "</a>" +
                                    "<button id='view_card' class='btn btn-sm btn-danger' name='"+data+"'><i class='glyphicon glyphicon-book'></i></button>" +
                                    "<button id='add_to_clipboard' class='btn btn-sm btn-success' name='"+data+"'><i class='glyphicon glyphicon-plus'></i></button></div>"
                            }
                        }
                    ],
                    "buttons": [
                        "copy"
                    ],
                    "searching": false,
                    "lengthChange": false,
                    "pageLength": 10
                });

    $("#download_json").attr("href", url);
}