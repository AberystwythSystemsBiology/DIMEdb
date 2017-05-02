$(document).ready(function () {
    $("input[type=radio][name=ionisation]").change(function () {
        var iad = {
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
            "positive": [["[M+H]1+", 1], ["[M+K]1+", 1], ["[2M+H]1-", 0]],
            "neutral": [["[M]", 1]]
        };

        $("#ionisation-adducts").empty();

        for (x in iad[this.value]) {
            var adduct = (iad[this.value][x]);
            if (adduct[1] == 1) {
                $("#ionisation-adducts").append("<option selected>" + adduct[0] + "</option>")
            }
            else {
                $("#ionisation-adducts").append("<option>" + adduct[0] + "</option>")
            }
        }
    });

    $('#advanced_search_results').on('click', "#add_to_clipboard", function() {
        clipboard_hook($(this).attr("name"));
    });

    $('#advanced_search_results').on('click', "#view_card", function() {
        generate_card($(this).attr("name"))
    });

    $("#search_button").click(function () {
        var masses = $("#masses").tagsinput("items");
        var ionisation = $("input[type=radio][name=ionisation]:checked").val();

        $("#advanced_search_results").fadeOut("slow");
        $("#advanced_search_results").empty();
        masses.forEach(function (mass, index) {
            api_url = generate_api_url(mass, ionisation);
            populate_results(mass, ionisation, api_url);
        });
    });
});

function calculate_ppm(mass, ppm_tolerance) {
    var tolerance = (mass * ppm_tolerance) * 0.0001;
    return [mass - tolerance, mass + tolerance]
}

function generate_api_url(mass, ionisation) {
    var ppm_tolerance = $("#ppm_tolerance").val();
    var masses = calculate_ppm(parseFloat(mass), ppm_tolerance);

    var adducts = [];

    $("#ionisation-adducts option:selected").each(function () {
        var $this = $(this);
        if ($this.length) {
            var selText = $this.text();
            adducts.push(selText)
        }
    });


    var base_api = 'api/metabolites/?where={"adducts.'+ionisation+'" : {"$elemMatch" : {"type" : {"$in" : '+JSON.stringify(adducts)+'},';

    base_api += '"accurate_mass" : {"$gte" :'+masses[0]+' , "$lte" : '+masses[1]+'}}}}';

    var projection = '&projection={"name" : 1, "molecular_formula" : 1, "adducts.'+ionisation+'.$":1}&max_results=1000';

    /*
    var origins = $("input[name='origins']:checked").map(function(){
            return $(this).val();
        }).toArray();

    if (origins.length > 0) {
        base_api += "&origins__contains=";
        for (i in origins) {
                if (i == origins.length-1) {
                    base_api += origins[i];
                }
                else {
                    base_api += origins[i]+",";
                }
        }
    }
    return base_api;
    */

    return base_api+projection
}

function generate_results_table(mass) {
    var mass_table = "<table id='search_results_"+mass.replace(".", "_")+"' class='table table-striped table-responsive display'><thead>"+
                "<tr></tr></thead></table>";
    return mass_table;
}


function populate_results(mass, ionisation, api_url) {
    $("#advanced_search_results").append("<h3> "+mass+" m/z</h3>");
    $("#advanced_search_results").append(generate_results_table(mass));
    generate_search_results(mass, ionisation, api_url);
    $('#advanced_search_results').fadeIn("slow");

}

function generate_search_results(mass, ionisation, api_url) {
    console.log(api_url);
     $('#search_results_'+mass.replace(".", "_")+'').DataTable({
         "destroy" : true,
         "searching": false,
         "lengthChange": false,
         "pageLength": 5,
         "ajax" : {
            "url" : encodeURI(getBaseURL()+api_url),
             "dataSrc" : "_items"
         },
         "columns" : [
             {
                            "title": "Metabolite Name",
                            "width": "40%",
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
     });
}
