$(document).ready(function () {
    var current_url = window.location.href;
    $('#search_results').hide();
    $("#search_button").click(function (e) {

        var query = $("#search_text").val();
        if (query != "") {
            var search_type = $("#s_selector_button").text();
            if (search_type == "Name") {
                var url = current_url + "api/metabolites/?name__icontains=" + String(query);
                $('#search_results').DataTable({
                    "destroy": true,
                    "ajax": url,
                    "columns": [
                        {
                            "title": "Name",
                            "width": "70%",
                            "data": "name",
                            "render": function (data, type, row) {
                                return data

                            }
                        },
                        {
                            "title": "Molecular Formula",
                            "width": "10%",
                            "data": "molecular_formula",
                            "render": function (data, type, row) {
                                return data.replace(/([0-9]+)/g, '<sub>$1</sub>');
                            }
                        },
                        {
                            "title": "Accurate Mass (m/z)",
                            "data": "accurate_mass",
                            "width": "10%",
                            "render": function (data, type, row) {
                                return data.toFixed(4);
                            }
                        },
                        {
                            "title": "Actions",
                            "data": "id",
                            "width": "10%",
                            "render": function (data, type, row) {
                                var view_url = current_url + "view/" + data;
                                return "<a href='" + view_url + "' target='_blank'><button class='btn btn-sm btn-primary'>View</button></a>"
                            }
                        }
                    ],
                    "searching": false,
                    //"bSort" : false,
                    "lengthChange": false,
                    "pageLength": 10
                });
            }
            else if (search_type == "Adduct") {
                var ion = $('input[name="ioinisation"]:checked').val();
                if (jQuery.isNumeric(query)) {
                    var url = current_url + "api/adducts/?adducts__pi_ppm=" + String(ion) + "," + String(query) + ",20";
                    $('#search_results').DataTable({
                        "destroy": true,
                        "ajax": url,
                        "columns": [
                            {
                                "title": "Name",
                                "width": "50%",
                                "data": "name",
                                "render": function (data, type, row) {
                                    return data

                                }
                            },
                            {
                                "title": "Molecular Formula",
                                "width": "10%",
                                "data": "molecular_formula",
                                "render": function (data, type, row) {
                                    return data.replace(/([0-9]+)/g, '<sub>$1</sub>');
                                }
                            },
                            {
                                "title": "Accurate Mass (m/z)",
                                "data": "adducts",
                                "width": "10%",
                                "render": function (data, type, row) {
                                    return data["neutral"]["peaks"][0]["accurate_mass"].toFixed(3);
                                }
                            },
                            {
                                "title": "Adduct",
                                "data": "adducts",
                                "width": "10%",
                                "render": function (data, type, row) {
                                    return data[ion]["peaks"][0]["type"];
                                }
                            },
                            {
                                "title": "Adduct (m/z)",
                                "data": "adducts",
                                "width": "10%",
                                "render": function (data, type, row) {
                                    return data[ion]["peaks"][0]["accurate_mass"].toFixed(3);
                                }
                            },
                            {
                                "title": "Difference (m/z)",
                                "data": "adducts",
                                "width": "5%",
                                "render": function (data, type, row) {
                                    return Math.abs(data[ion]["peaks"][0]["accurate_mass"] - parseFloat(query)).toFixed(3);
                                }
                            },
                            {
                                "title": "Actions",
                                "data": "id",
                                "width": "5%",
                                "render": function (data, type, row) {
                                    var view_url = current_url + "view/" + data;
                                    return "<a href='" + view_url + "' target='_blank'><button class='btn btn-sm btn-primary'>View</button></a>"
                                }
                            }
                        ],
                        "searching": false,
                        //"bSort" : false,
                        "lengthChange": false,
                        "pageLength": 10
                    });
                }
            }
            else if (search_type == "Accurate Mass") {
                if (jQuery.isNumeric(query)) {
                var url = current_url + "api/adducts/?adducts__pi_ppm=neutral," + String(query) + ",20";
                $('#search_results').DataTable({
                        "destroy": true,
                        "ajax": url,
                        "columns": [
                            {
                                "title": "Name",
                                "width": "50%",
                                "data": "name",
                                "render": function (data, type, row) {
                                    return data

                                }
                            },
                            {
                                "title": "Molecular Formula",
                                "width": "10%",
                                "data": "molecular_formula",
                                "render": function (data, type, row) {
                                    return data.replace(/([0-9]+)/g, '<sub>$1</sub>');
                                }
                            },
                            {
                                "title": "Accurate Mass (m/z)",
                                "data": "adducts",
                                "width": "10%",
                                "render": function (data, type, row) {
                                    return data["neutral"]["peaks"][0]["accurate_mass"].toFixed(3);
                                }
                            },
                            {
                                "title": "Difference (m/z)",
                                "data": "adducts",
                                "width": "5%",
                                "render": function (data, type, row) {
                                    return Math.abs(data["neutral"]["peaks"][0]["accurate_mass"] - parseFloat(query)).toFixed(3);
                                }
                            },
                            {
                                "title": "Actions",
                                "data": "id",
                                "width": "5%",
                                "render": function (data, type, row) {
                                    var view_url = current_url + "view/" + data;
                                    return "<a href='" + view_url + "' target='_blank'><button class='btn btn-sm btn-primary'>View</button></a>"
                                }
                            }
                        ],
                        "searching": false,
                        //"bSort" : false,
                        "lengthChange": false,
                        "pageLength": 10
                    });
                }


            }
            $("#search_results").fadeIn("slow");
        }
        else {
        }
    });

    $('#search_text').keypress(function (e) {
        if (e.which == 13) {//Enter key pressed
            $('#search_button').click();//Trigger search button click event
        }
    });
});

$("#s_selector li").click(function () {
    var v = $(this).text();
    $("#s_selector_button").text(v);
    if (v == "Name") {
        $("#hp_ionisation_selector").fadeOut("slow");
    }

    if (v == "Accurate Mass") {
        $("#hp_ionisation_selector").fadeOut("slow");
    }

    if (v == "Adduct") {
        $("#hp_ionisation_selector").fadeIn("slow");
    }

});