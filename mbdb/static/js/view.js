function render_view(base_url, id) {
    var url = base_url + "api/metabolite/?id__exact=" + id;
    var x_plot = [];
    var y_plot = [];

    $.getJSON(url, function (data) {
        var metabolite = data["data"][0];
        // Change the document title.
        $(document).prop('title', metabolite["name"] + " : MetaBolomics DataBase");
        $("#met_name").html(metabolite["name"]);
        $("#get_json").attr("href", base_url+"api/metabolite/?id__exact="+metabolite["id"]);
        $("#molecular_formula").html(metabolite["molecular_formula"].replace(/([0-9]+)/g, '<sub>$1</sub>'));
        $("#molecular_formula_search").attr("href", base_url+"api/metabolites/?molecular_formula__exact="+metabolite["molecular_formula"])
        $("#accurate_mass").html(metabolite["accurate_mass"]);
        $("#neutral_mass").html(metabolite["adduct_weights"]["neutral"]);
        // $("#origins").html(metabolite["origins"]);
        $("#smiles").html(metabolite["smiles"]);

        for (o in metabolite["origins"]) {
            $("#origins").append("<p>"+metabolite["origins"][o]+"</p>")
        }

        for (result in metabolite["adduct_weights"]["negative"]["peaks"]) {
            var peak = metabolite["adduct_weights"]["negative"]["peaks"][0]
            $("#negative_adduct").append("<li class='list-group-item'><b>"+peak[0]+":</b> "+peak[1].toFixed(4)+"</li>")
        }

        $.ajax({
            url: base_url+"gen_structure/"+metabolite["id"],
            success: function (result) {
                $("#structure").attr("src", "data:image/png;base64,"+result);
            }
        });
        for (i in metabolite["isotopic_distributions"]) {
            var spectra = metabolite["isotopic_distributions"][i];
            var mz_spectra = (spectra[0] + metabolite["adduct_weights"]["neutral"]);
            x_plot.push(mz_spectra);
            y_plot.push(spectra[1]);
            $("#distribution_table tbody").append("<tr><td>"+mz_spectra.toFixed(4)+"</td><td>"+spectra[1].toFixed(2)+"</td></tr>")
        }


        var isotopic_data = [{
            x: x_plot,
            y: y_plot,
            type: 'bar'
        }];

        console.log();
        console.log(y_plot);

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

        Plotly.newPlot("distribution_chart", isotopic_data, layout)
        $("#container").show()
    });
}