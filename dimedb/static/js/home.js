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
    var search_query = '?where={"name" : {"$regex": ".*'+query+'.*"}}&projection={"name" : 1, "accurate_mass" : 1, "molecular_formula" : 1}&max_results=1000';

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
                "title": "Molecular Formula",
                "className": "dt-center",
                "width" : "15%",
                "data": "molecular_formula",
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