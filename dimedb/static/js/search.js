function clear_results() {
    $("#search_results").fadeOut("slow");
    $("#search_results").empty();
}


function generate_table(mass, ionisation, api_url) {

    function insert_table(mass) {
        var mass_table = "<table id='results_" + mass.replace(".", "_") + "' class='table table-striped table-responsive display' width='100%'><thead>" +
        "<tr></tr></thead></table>";
        return mass_table
    }

    function generate_datatable(mass, ionisation, api_url) {
        $('#results_' + mass.replace(".", "_") + '').DataTable({
        "destroy": true,
        "searching": false,
        "lengthChange": false,
        "pageLength": 5,
        "ajax": {
            "url": encodeURI(api_url),
            "dataSrc": "_items"
        },
        "columns": [
            {
                "title": "Metabolite Name",
                "width": "40%",
                "data": "name",
                "render": function (data, type, row) {
                    return "<a href='"+getBaseURL()+"view/"+row._id+"'>"+row.name+"</a>"
                }
            },
            {
                "title": "Molecular Formula",
                "width": "10%",
                "className": "dt-center",
                "data": "molecular_formula",
                "render": function (data, type, row) {
                    return '<p style="text-align:center">' + data.replace(/([0-9]+)/g, '<sub>$1</sub>') + '</p>';
                }
            },
            {
                "title": "Mass (m/z)",
                "data": "adducts",
                "className": "dt-center",
                "width": "10%",
                "render": function (data, type, row) {
                    return data[ionisation][0]["accurate_mass"];
                }
            },
            {
                "title": "Adduct",
                "data": "adducts",
                "className": "dt-center",
                "width": "10%",
                "render": function (data, type, row) {
                    return data[ionisation][0]["type"];
                }
            },
            {
                "title": "Difference",
                "data": "adducts",
                "className": "dt-center",
                "width": "10%",
                "render": function (data, type, row) {
                    return (data[ionisation][0]["accurate_mass"] - mass).toFixed(4);
                }
            },
            {
                "title": "Actions",
                "data": "_id",
                "className": "dt-center",
                "width": "15%",
                "render": function (data, type, row) {

                    var row_info = "<div id='tester' class='btn-toolbar text-center'>";
                    if (ionisation == "neutral") {
                        row_info += "<button id='view_card' class='btn btn-sm btn-danger' name='" + data + "," + row["adducts"]["neutral"][0]["type"] +  "'>"
                    }
                    else if (ionisation == "negative") {
                        row_info += "<button id='view_card' class='btn btn-sm btn-danger' name='" + data + "," + row["adducts"]["negative"][0]["type"] +  "'>"
                    }

                    else {
                        row_info += "<button id='view_card' class='btn btn-sm btn-danger' name='" + data + "," + row["adducts"]["positive"][0]["type"] +  "'>"
                    }

                    row_info += "<i class='glyphicon glyphicon-book'></i></button>";
                    row_info += "<button id='add_to_clipboard' class='btn btn-sm btn-success' name='" + data + "'><i class='glyphicon glyphicon-plus'></i></button></div>";

                    return row_info

                }
            }
        ],
});
    }

    $("#search_results").append("<h3> " + mass + " m/z</h3>");
    $("#search_results").append(insert_table(mass));
    generate_datatable(mass, ionisation, api_url);
}

function generate_api_url(mass, ionisation) {

    function generate_ppm(mass) {
        var ppm_tolerance = $("#ppm_tolerance").val();
        return (mass / ((1/ppm_tolerance)*10000))
    }

    function generate_adducts() {
        var adducts = [];
        $('#ionisation_adducts').find(":selected").each(function() {
            adducts.push($(this).val());
        });
        return JSON.stringify(adducts)
    }

    function generate_origins() {
        var selected_origins = [];
        $('#origins input:checked').each(function() {
            selected_origins.push($(this).attr("value"));
        });
        return selected_origins
    }

    function generate_biofluids() {
        var selected_biofluids = [];
        $('#biofluid_location input:checked').each(function() {
            selected_biofluids.push($(this).attr("value"));
        });
        return selected_biofluids
    }

    function generate_tissues() {
        var selected_tissues = [];
        $('#tissue_location input:checked').each(function() {
            selected_tissues.push($(this).attr("value"));
        });
        return selected_tissues
    }

    var tolerance = generate_ppm(mass);

    var api_url = getBaseURL()+"api/metabolites/?where={";




    api_url += '"adducts.' + ionisation + '" : {"$elemMatch" : {"type" : {"$in" : ' + generate_adducts() + '},';


    var lte = parseFloat(mass) + parseFloat(tolerance);
    var gte = parseFloat(mass) - parseFloat(tolerance);


    api_url += '"accurate_mass" : {"$lte" : ' + lte +  ', "$gte" : ' + gte +'}}}';

    var origins = generate_origins();


    if (origins.length > 0) {
        api_url += ',"origins" : { "$in" : ' + JSON.stringify(origins) + '}';
    }

    var biofluids = generate_biofluids();


    if (biofluids.length > 0) {
        api_url += ',"biofluid_location" : { "$in" : ' + JSON.stringify(biofluids) + '}';
    }

    var tissues = generate_tissues();

    if (tissues.length > 0) {
        api_url += ',"tissue_locations" : { "$in" : ' + JSON.stringify(tissues) + '}';

    }


    api_url += '}&projection={"name" : 1, "molecular_formula" : 1, "adducts.' + ionisation + '.$":1}&max_results=1000}';
    return api_url;
}

$("#search_button").click(function () {
    clear_results();
    var masses = $("#masses").tagsinput("items");
        var ionisation = $("input[type=radio][name=ionisation]:checked").val();
        masses.forEach(function(mass, index) {
            var mass = mass.trim();
            var api_url = generate_api_url(mass, ionisation);
            generate_table(mass, ionisation, api_url);
        });

    $("#search_results").fadeIn("slow");
});

$(document).ready(function () {
    function fill_adducts_list(ionisation) {
        var adducts_dict = {
            "negative": [
                ["[M-H]1-", 1],
                ["[M+Cl]1-", 1],
                ["[M+Br]1-", 1],
                ["[M+Na-2H]1-", 0],
                ["[M+K-2H]1-", 0],
                ["[3M-H]1-", 0],
                ["[2M+Hac-H]1-", 0],
                ["[2M+FA-H]1-", 0],
                ["[2M-H]1-", 0],
                ["[M+TFA-H]1-", 0],
                ["[M+Hac-H]1-", 0],
                ["[M+FA-H]1-", 0],
                ["[M-2H]2-", 0],
                ["[M-3H]3-", 0],
                ["[2M+Na-2H]1-", 0],
                ["[M1-.]1-", 0]
            ],
            "positive": [
                ["[M+H]1+", 1],
                ["[M+K]1+", 1],
                ["[M+Na]1+", 1],
                ["[2M+H]1-", 0]
            ],
            "neutral": [
                ["[M]", 1]
            ]
        };

        $("#ionisation_adducts").empty();

        adducts_dict[ionisation].forEach(function(adduct, index) {
            if (adduct[1] == 1) {
                $("#ionisation_adducts").append($("<option selected></option>").text(adduct[0]).attr("value", adduct[0]));
            }
            else {
                $("#ionisation_adducts").append($("<option></option>").text(adduct[0]).attr("value", adduct[0]));
            }
        });


    }

    $("input[type=radio][name=ionisation]").change(function () {
        fill_adducts_list(this.value)
    });

});

