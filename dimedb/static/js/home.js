function naming_is_hard(query) {
    generate_homepage_results(query);
    $("#home-page-welcome-search").fadeOut();
    $("#results_search_text").val(query);
    $("#home-page-search-results").delay(460).fadeIn();
}


$(document).ready(function () {

    $(".navbar img").hide();

    $('#home-page-search-results').hide();

    $("#search_button").click(function (e) {
        var query = $("#search_text").val();
        if (query != "") {
            naming_is_hard(query);
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

    $("#results_search_button").click(function (e) {
        var query = $("#results_search_text").val();
        if (query != "") {
            naming_is_hard(query);
        }
        else {
            alert("You don't seem to have entered anything?");
        }
    });

    $('#results_search_text').keypress(function (e) {
        if (e.which == 13) {
            $('#results_search_button').click();
        }
    });

    $("#search_button").click(function (e) {
        var query = $("#search_text").val();
        if (query != "") {
            $("#search_results").fadeOut();
            generate_homepage_results(query);
            $("#search_results").delay(600).fadeIn();

        }
        else {
            alert("You don't seem to have entered anything?");
        }
    });

    $('#search_text_results').keypress(function (e) {
        if (e.which == 13) {
            $('#search_button').click();
        }
    });

    $('#home-page-search-results').on('click', "#view_card", function() {
        generate_card($(this).attr("name"))
    });

    $('#home-page-search-results').on('click', "#add_to_clipboard", function() {
        clipboard_hook($(this).attr("name"));
    });

});

function generate_homepage_results(query) {
    var search_query = '?where={"name" : {"$regex": ".*'+query+'.*"}}&projection={"name" : 1, "accurate_mass" : 1, "chemical_formula" : 1}&max_results=1000';

    console.log(search_query);

    $('#search_results').DataTable({
        "destroy" : true,
        "searching": false,
        "lengthChange": false,
        "fixedColumns": true,
        "pageLength": 30,
        "ajax" : {
            "url" : encodeURI(getBaseURL()+"api/metabolites/"+search_query),
            "dataSrc" : "_items"
        },
        "columns" : [
            {
                "title": "Metabolite Name",
                "data": "name",
                "width" : "75%",
                "render": function (data, type, row) {
                    return "<a href='"+getBaseURL()+"view/"+row._id+"'>"+row.name+"</a>"
                }
            },
            {
                "title": "Chemical Formula",
                "className": "dt-center",
                "width" : "15%",
                "data": "chemical_formula",
                "render": function (data, type, row) {
                    return data.replace(/([0-9]+)/g, '<sub>$1</sub>');
                }
            },
            {
                "title": "Mass (m/z)",
                "className": "dt-center",
                "width" : "7.5%",
                "data": "accurate_mass",
                "render": function (data, type, row) {
                    return data.toFixed(4);
                }
            },
            {
                "title": "Actions",
                "data": "_id",
                "width" : "10%",
                "className": "dt-center",
                "render": function (data, type, row) {
                    var view_url = getBaseURL() + "view/" + data;
                    return "<div id='tester' class='btn-toolbar text-center'>" +
                        "<button id='view_card' class='btn btn-sm btn-danger' name='"+data+"'><i class='glyphicon glyphicon-book'></i></button>" +
                        "<button id='add_to_clipboard' class='btn btn-sm btn-success' name='"+data+"'><i class='glyphicon glyphicon-plus'></i></button></div>"
                }
            }
        ],
    });
}

function generate_card(id) {


    function fill_origins(origins_array) {
        if (origins_array.length != 0) {
            var origins = "";
            for (i = 0; i < origins_array.length; i++) {
                origins += origins_array[i] + "; ";
            }
        }

        else {
            var origins = "Not available";
        }

        return origins;
    }

    function fill_biofluids(biofluids_array) {
        if (biofluids_array != null) {
            var biofluids = "";
            for (i = 0; i < biofluids_array.length; i++) {
                biofluids += biofluids_array[i] + "; ";
            }
        }

        else {
            var biofluids = "Not available";
        }

        return biofluids;
    }

    function fill_tissue_locations(tissue_locations_array) {
        if (tissue_locations_array != null) {
            var tissues = "";
            for (i = 0; i < tissue_locations_array.length; i++) {
                tissues += tissue_locations_array[i] + "; ";
            }
        }
        else {
            var tissues = "Not available";
        }

        return tissues
    }

    function fill_adduct_panels(adducts) {
        var ionisations = ["negative", "positive"];

        for (i in ionisations) {
            var ionisation = ionisations[i];
            for (j in adducts[ionisation]) {
                var adduct  = adducts[ionisation][j];
                var name = ionisation + "," + adduct["type"];
                $("#card_"+ionisation+"_adduct").append(
                    "<li class='list-group-item'><b>"+ adduct["type"] + ":</b> "+adduct["accurate_mass"].toFixed(4)+ "</li>"
                )
            }
        }

    }

    var metabolite = get_metabolite(id);

    console.log(metabolite);

    $("#card_metabolite_name").html(metabolite["name"]);
    $("#card_molecular_formula").html(metabolite["chemical_formula"].replace(/([0-9]+)/g, '<sub>$1</sub>'));
    $("#card_accurate_mass").html(metabolite["accurate_mass"].toFixed(6));
    $("#card_neutral_mass").html(metabolite["adducts"]["neutral"][0]["accurate_mass"].toFixed(6));
    $("#card_origins").html(fill_origins(metabolite["origins"]));
    $("#card_biofluids").html(fill_biofluids(metabolite["biofluid_locations"]));
    $("#card_tissues").html(fill_tissue_locations(metabolite["tissue_locations"]));
    fill_adduct_panels(metabolite["adducts"]);
    $("#metabolite_card").modal("toggle");
}