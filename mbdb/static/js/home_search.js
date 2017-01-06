$(document).ready(function () {
    var current_url = window.location.href;
    $("#search_results").hide();
    $("#search_button").click(function (e) {
        query = $("#search_text").val();
        if (query != "") {
            var url = current_url + "api/metabolites/?name__contains="+String(query);
            $("#search_results").show();

            $('#search_results').DataTable({
                "ajax": url,
                "columns": [
                    {
                        "data": "name",
                        "render": function (data, type, row) {
                            return data

                        }
                    },
                    {
                        "data": "molecular_formula",
                        "render": function (data, type, row) {
                            return data.replace(/([0-9]+)/g, '<sub>$1</sub>');
                        }
                    },
                    {
                        "data": "accurate_mass",
                        "render": function (data, type, row) {
                            return data.toFixed(4);
                        }
                    },
                    {
                        "data" : "id",
                        "render" : function (data, type, row) {
                            return "<button class='btn btn-sm btn-primary disabled'>View</button>"
                        }
                    }
                ],
                "searching": false,
                "bSort" : false,
                "lengthChange": false,
                "pageLength": 10,
                "bRetrieve" : true
            });


        }
        else {
            alert("No query entered?")
        }
    });
});