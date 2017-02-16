$(document).ready(function () {
    $("input[type=radio][name=ionisation]").change(function() {

        var iad = {
            "negative" : [["M-H", 1], ["M+Cl", 1], ["M+Br", 0]],
            "positive" : [["M+H", 1], ["M+K", 1], ["2M+H", 0]],
            "neutral" : [["M", 1]]
        }

        $("#ionisation-adducts").empty()

        for(x in iad[this.value]) {
            var adduct = (iad[this.value][x]);
            if (adduct[1] == 1) {
                $("#ionisation-adducts").append("<option selected>"+adduct[0]+"</option>")
            }
            else {
                $("#ionisation-adducts").append("<option>"+adduct[0]+"</option>")
            }
        }

    });

    $("#search_button").click(function() {
        var ionisation = $("input[type=radio][name=ionisation]:checked").val();
        var current_url = window.location.protocol + "//" + window.location.host

        var base_api = "/api/adducts/?adducts__pi_ppm=";
        base_api += ionisation + ",";
        base_api += $("#mass").val() +",";
        base_api += $("#ppm_tolerance").val();

        var origins = $("input[name='origins']:checked").map(function(){
            return $(this).val();
        }).toArray();

        var biofluids = $("input[name='biofluids']:checked").map(function(){
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

        if (biofluids.length > 0) {
            base_api += "&biofluids__contains=";
            for (i in biofluids) {
                if (i == biofluids.length-1) {
                    base_api += biofluids[i];
                }
                else {
                    base_api += biofluids[i]+",";
                }
            }
        }

        render_search_results(base_api);
    });

});

function render_search_results(url) {
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
                    "pageLength": 10
                });

}