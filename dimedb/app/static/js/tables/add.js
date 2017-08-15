function get_usertables() {
    var json = (function () {
        var json = null;
        $.ajax({
            'async': false,
            'global': false,
            'url': getBaseURL()+"tables/mytables/api/get_mytables",
            'dataType': "json",
            'success': function (data) {
                json = data;
            }
        });
        return json;
    })();
    return json["data"]
}

function populate_select_area(tables) {

    if (tables.length > 0) {
        $("#no_table_warning").fadeOut(0);
        $("#add_to_table_container").fadeIn(200);
        for (index in tables) {
            var option = new Option("DdbT" + tables[index]["id"] + ": " + tables[index]["Title"], tables[index]["id"])
            $("#table_selector").append(option);
        }
    }
    else {
        $("#add_to_table_container").fadeOut(200);

    }

}

$(document).ready(function () {
    $("#success").fadeOut(0);
    var tables = get_usertables();
    populate_select_area(tables);

    $("#add_to_table").click(function () {
        $("#add_table_modal").modal().toggle()
    });
});

$("#add_submission").click(function () {
   var payload = {
       "Table ID" : $("#table_selector").find(":selected").val(),
       "InChI Key" : $("#inchi_key").text(),
       "Comments" : $("#comments").val()
   };

   $.ajax({
       "type" : "POST",
       "url" : getBaseURL()+"tables/mytables/api/add_metabolite",
       "data" : payload,
       "success" : function () {
           $("#add_to_table_container").fadeOut(200);
           $("#success").fadeIn(500);
       }
   })
});