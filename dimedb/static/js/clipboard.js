

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

function change_count() {
}

function inArr(val) {
    val = JSON.stringify(val);
    return metabolite_array.filter(function (el) {
            return JSON.stringify(el) == val
        }
    ).length;
}

function metabolite_to_clipboard(id) {
    var clipboard_array = [id, metabolite["name"], metabolite["sources"]["kegg_id"]];

    if (inArr(clipboard_array) != 0) {
        alert("Metabolite already in clipboard");
    }
    else {
        metabolite_array.push(clipboard_array);
        populate_clipboard(clipboard_array);
        localStorage.setItem("metabolite_array", JSON.stringify(metabolite_array));
        change_count();
    }
}

function populate_clipboard(m) {
    console.log(m);
}

$(document).ready(function () {
    var metabolite_array = load_metabolites();
    $("#mc_count").html(metabolite_array.length)
});

$('#metabolite_clipboard_panel').on('click', "#clear_clipboard", function () {
    $("#clipboard_list_group li").remove();
    localStorage.clear("metabolite_array");
});




