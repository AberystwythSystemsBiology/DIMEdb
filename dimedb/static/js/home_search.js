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

    $("#search_button").click(function (e) {
        var query = $("#search_text").val();
        if (query != "") {
            $("#search_results").fadeOut();
            render_search_results("search_results", "api/metabolites/?name__icontains="+query, 10);
            $("#search_results").delay(460).fadeIn();

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