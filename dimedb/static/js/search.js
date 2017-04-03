$(document).ready(function () {


    $("input[type=radio][name=ionisation]").change(function () {
        var iad = {
            "negative": [["M-H", 1], ["M+Cl", 1], ["M+Br", 0]],
            "positive": [["M+H", 1], ["M+K", 1], ["2M+H", 0]],
            "neutral": [["M", 1]]
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

    $("#search_button").click(function () {
        var masses = $("#masses").tagsinput("items");
        $("#advanced_search_results").fadeOut("slow");
        $("#advanced_search_results").empty();
        masses.forEach(function (mass, index) {
            api_url = generate_api_url(mass);
            populate_results(mass, getBaseURL()+api_url);
        });
    });
});


function generate_api_url(mass) {
    var ionisation = $("input[type=radio][name=ionisation]:checked").val();
    var ppm_tolerance = $("#ppm_tolerance").val();

    var base_api = "api/adducts/?adducts__pi_ppm=";
    base_api += ionisation+ "," + mass + "," + ppm_tolerance;
    return base_api;
}

function generate_results_table(mass) {
    var mass_table = "<table id='search_results_"+mass.replace(".", "_")+"' class='table table-striped table-responsive display'><thead>"+
                "<tr></tr></thead></table>";
    return mass_table;
}


function populate_results(mass, api_url) {
    $("#advanced_search_results").append("<h3> "+mass+" m/z</h3>");
    $("#advanced_search_results").append(generate_results_table(mass));
    render_search_results(mass, api_url);
    $('#advanced_search_results').fadeIn("slow");

}


function render_search_results(mass, api_url) {
    console.log(api_url);
    $('#search_results_'+mass.replace(".", "_")).delay(100).DataTable({
                    "destroy": true,
                    "ajax": api_url,
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
                            "width": "5%",
                            "render": function (data, type, row) {
                                var view_url = window.location.protocol + "//" + window.location.host + "/" + "view/" + data;
                                return "<a href='" + view_url + "' target='_blank'><button class='btn btn-sm btn-primary' id='view_button'>View</button></a>"
                            }
                        }
                    ],
                    "searching": false,
                    //"bSort" : false,
                    "lengthChange": false,
                    "pageLength": 5
                });

}

function getBaseURL () {
   return location.protocol + "//" + location.hostname +
      (location.port && ":" + location.port) + "/";
}
