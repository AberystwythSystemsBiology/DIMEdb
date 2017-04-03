var base_url = window.location.protocol + "//" + window.location.host + "/" + window.location.pathname.split('/')[1];
var metabolite_array = load_metabolites();

$(document).ready(function () {
    change_count();

    $('#metabolite_clipboard_panel').on('click', "#clear_clipboard", function () {
        $("#clipboard_list_group li").remove();
        localStorage.clear("metabolite_array");
    });

    $("#metabolite_clipboard_panel").on("click", "#clipboard_to_kegg", function () {
        kegg_map_url = "http://www.genome.jp/kegg-bin/show_pathway?map01100+";
        metabolite_array.forEach(function (metabolite, index) {
            if (metabolite[2] != null) {
                $("#pathway_mapper_h_buttons").append(
                    "<a href='http://www.genome.jp/dbget-bin/www_bget?compound+" + metabolite[2] + "' target='_blank'>" +
                    "<button class='btn btn-default'>" + metabolite[2] + "</button></a>"
                );
                kegg_map_url += metabolite[2] + "+";
            }
        });
        kegg_map_url = kegg_map_url.slice(0, -1);
        $("#mapped_pathway_iframe").attr("src", kegg_map_url);
        $("#view_mapped_pathway").modal("toggle");
    });

});

function change_count() {
    $("#clipboard_count").html(metabolite_array.length)
}

function inArr(val) {
    val = JSON.stringify(val);
    return metabolite_array.filter(function (el) {
            return JSON.stringify(el) == val
        }
    ).length;
}

function clipboard_hook(id) {
    get_metabolite(id, function (metabolite) {
        var m = [id, metabolite["name"], metabolite["sources"]["kegg_id"]];

        if (inArr(m) != 0) {
            alert("Metabolite already in clipboard");
        }
        else {
            metabolite_array.push(m);
            populate_clipboard(m);
            localStorage.setItem("metabolite_array", JSON.stringify(metabolite_array));
            change_count();
        }
    });
}

function populate_clipboard(m) {
    $("#clipboard_list_group").append(
        "<li id='" + m[0] + "' class='list-group-item'>" + m[1] + "<a href='" + base_url + "view/" + m[0] + "'>" +
        "<button class='btn btn-success btn-sm pull-right'>" +
        "<i class='glyphicon glyphicon-link'></i> View" +
        "</button></a><div class='clearfix'></div></li>"
    );
}

function get_metabolite(id, callback) {
    $.getJSON(base_url + "api/metabolite/?id=" + id, function (data) {
        callback(data["data"][0]);
    });
}

function load_metabolites() {
    var metabolite_array = [];
    var stored_metabolites = JSON.parse(localStorage.getItem("metabolite_array"));
    if (stored_metabolites != null) {
        metabolite_array = stored_metabolites;
        metabolite_array.forEach(function (metabolite, index) {
            populate_clipboard(metabolite);
        });
    }
    return metabolite_array
}