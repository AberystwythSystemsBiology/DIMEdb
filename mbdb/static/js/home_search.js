$(document).ready(function () {
    var current_url = window.location.href;
    $('#search_results').hide()
    $("#search_button").click(function (e) {
        var query = $("#search_text").val();
        if (query != "") {
            var url = current_url + "api/metabolites/?name__icontains="+String(query);
            $('#search_results').DataTable({
                "destroy" : true,
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
                            var view_url = current_url+"view/"+data;
                            return "<a href='"+view_url+"' target='_blank'><button class='btn btn-sm btn-primary'>View</button></a>"
                        }
                    }
                ],
                "searching": false,
                //"bSort" : false,
                "lengthChange": false,
                "pageLength": 10
            });
            $('#search_results').show()
        }
        else {
        }
    });

    $('#search_text').keypress(function(e){
        if(e.which == 13){//Enter key pressed
            $('#search_button').click();//Trigger search button click event
        }
    });
});