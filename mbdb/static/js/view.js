function render_view(base_url, id) {
    var url = base_url + "api/metabolite/?id__exact=" + id;
    var x_plot = [];
    var y_plot = [];

    $.getJSON(url, function (data) {
        // Obtain metabolite dictonary
        var metabolite = data["data"][0];

        // Change the document title.
        $(document).prop('title', metabolite["name"] + " : MetaBolomics DataBase");

        // Fill up the Information Table.
        $("#met_name").html(metabolite["name"]);
        $("#get_json").attr("href", base_url+"api/metabolite/?id__exact="+metabolite["id"]);
        $("#molecular_formula").html(metabolite["molecular_formula"].replace(/([0-9]+)/g, '<sub>$1</sub>'));
        $("#molecular_formula_search").attr("href", base_url+"api/metabolites/?molecular_formula__exact="+metabolite["molecular_formula"])
        $("#accurate_mass").html(metabolite["accurate_mass"]);
        $("#neutral_mass").html(metabolite["adducts"]["neutral"]["peaks"][0]["accurate_mass"]);
        $("#smiles").html(metabolite["smiles"]);

        for (o in metabolite["origins"]) {
            $("#origins").append("<p>"+metabolite["origins"][o]+"</p>");
        }


        $.ajax({
            url: base_url+"gen_structure/"+metabolite["id"],
            success: function (result) {
                $("#structure").attr("src", "data:image/png;base64,"+result);
            }
        });

        // Negative
        for (result in metabolite["adducts"]["negative"]["peaks"]) {
            var peak = metabolite["adducts"]["negative"]["peaks"][result];
            $("#negative_adduct").append("<li class='list-group-item'><b>"+peak["type"]+":</b> "+peak["accurate_mass"].toFixed(4)+"</li>");
        }

        for (result in metabolite["adducts"]["positive"]["peaks"]) {
            var peak = metabolite["adducts"]["positive"]["peaks"][result];
            $("#positive_adduct").append("<li class='list-group-item'><b>"+peak["type"]+":</b> "+peak["accurate_mass"].toFixed(4)+"</li>");
        }



        for (i in metabolite["adducts"]["neutral"]["peaks"][0]["isotopic_distribution"]) {
            var spectra = metabolite["adducts"]["neutral"]["peaks"][0]["isotopic_distribution"][i];
            var mz_spectra = (spectra[0]);
            x_plot.push(mz_spectra);
            y_plot.push(spectra[1]);
            $("#distribution_table tbody").append("<tr><td>"+mz_spectra.toFixed(4)+"</td><td>"+spectra[1].toFixed(2)+"</td></tr>")
        }


        var isotopic_data = [{
            x: x_plot,
            y: y_plot,
            type: 'bar'
        }];

        var layout = {
            xaxis: {
                title: 'Mass (mz)',
                showgrid: false,
                range: [Math.min.apply(Math, x_plot)-10, Math.max.apply(Math, x_plot)+10]
            },
            yaxis: {
                title: 'Intensity (%)',
                showline: false
            },
            margin: {
                l: 50,
                r: 50,
                b: 50,
                t: 50,
                pad: 4
            }
        };

        Plotly.newPlot("distribution_chart", isotopic_data, layout, {displayModeBar : false})
        $("#container").show()
    });
}