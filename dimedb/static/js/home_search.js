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
    })

    $('#search_text_results').keypress(function (e) {
        if (e.which == 13) {
            $('#search_text_button').click();
        }
    });


});


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
                            "width": "5%",
                            "render": function (data, type, row) {
                                var view_url = current_url + "view/" + data;
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