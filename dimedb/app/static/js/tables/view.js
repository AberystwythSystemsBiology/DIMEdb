$(document).ready(function () {
});


function generate_table(id) {

    $('#metabolite_table').DataTable({
            "destroy": true,
            "sort" : false,
            "pageLength": 10,
            "ajax": {
                "url": encodeURI(getBaseURL()+"tables/DdbT"+id+"/api/get_metabolites"),
                "dataSrc": "data"
            },
            "columns": [
                {
                    "title" : "Molecular Structure",
                    "width": "10%",
                    "render" : function(data,type,row) {
                        return "<img src='" + getBaseURL() + "view/structure/" + row.InChIKey + "' class='img-responsive img-circle'>"                    }
                },
                {
                    "title" : "Metabolite Name",
                    "width" : "80%",
                    "render" : function (data, type, row) {
                        var s = "<a href='" + getBaseURL() + "view/" + row.InChIKey + "' target='_blank'>" + row.Name + "</a>";
                        s += "<p>"+row.Comments+"</p>"
                        return s
                    }
                },

                {
                    "title" : "Molecular Formula",
                    "width" : "10%",
                    "className" : "dt-center",
                    "render" : function (data, type, row) {
                        return row["Molecular Formula"].replace(/([0-9]+)/g, '<sub>$1</sub>');
                    }
                }
            ],
        });
}
