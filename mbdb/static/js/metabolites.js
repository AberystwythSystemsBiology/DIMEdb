$(document).ready(function() {
    $('#metabolites').DataTable( {
        "ajax" : "http://127.0.0.1:5000/met/?_fields=id,name,molecular_formula,accurate_mass",
        "columns": [
            { "data" : "name",
            "render" : function (data, type, row) {
                return data

            }},
            { "data" : "molecular_formula",
                "render" : function(data, type, row) {
                    return data.replace(/([0-9]+)/g, '<sub>$1</sub>');
                }
            },
            { "data" : "accurate_mass",
                "render" : function (data, type, row) {
                    return data.toFixed(4);
                }
            }
        ],
        "searching" : false,
        "lengthChange" : false,
        "pageLength" : 50
    });
} );