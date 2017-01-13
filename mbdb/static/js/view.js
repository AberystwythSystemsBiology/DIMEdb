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
        $("#get_json").attr("href", base_url + "api/metabolite/?id__exact=" + metabolite["id"]);
        $("#molecular_formula").html(metabolite["molecular_formula"].replace(/([0-9]+)/g, '<sub>$1</sub>'));
        $("#molecular_formula_search").attr("href", base_url + "api/metabolites/?molecular_formula__exact=" + metabolite["molecular_formula"])
        $("#accurate_mass").html(metabolite["accurate_mass"]);
        $("#smiles").html(metabolite["smiles"]);

        for (o in metabolite["origins"]) {
            $("#origins").append("<p>" + metabolite["origins"][o] + "</p>");
        }


        $.ajax({
            url: base_url + "gen_structure/" + metabolite["id"],
            success: function (result) {
                $("#structure").attr("src", "data:image/png;base64," + result);
            }
        });

        for (result in metabolite["adducts"]) {
            if (result == "neutral") {
                $("#neutral_mass").html(metabolite["adducts"]["neutral"]["peaks"][0]["accurate_mass"]);
            }

            else {
                for (adduct_idx in metabolite["adducts"][result]["peaks"]) {
                    var adduct = metabolite["adducts"][result]["peaks"][adduct_idx];
                    $("#"+result+"_adduct").append("<li class='list-group-item'><b>"+adduct["type"]+":</b> "+adduct["accurate_mass"].toFixed(4)+"<button id='isotope' name='"+result+"_"+adduct_idx+"' class='btn btn-primary btn-sm pull-right'>Isotope</button><div class='clearfix'></div></li>");
                }
            }
        }

        for (i in metabolite["adducts"]["neutral"]["peaks"][0]["isotopic_distribution"]) {
            var spectra = metabolite["adducts"]["neutral"]["peaks"][0]["isotopic_distribution"][i];
            x_plot.push(spectra[0]);
            y_plot.push(spectra[1]);
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
                range: [Math.min.apply(Math, x_plot) - 10, Math.max.apply(Math, x_plot) + 10]
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

        Plotly.newPlot("distribution_chart", isotopic_data, layout, {displayModeBar: false});
        $("#container").show();

        // Prep table

        $('[id="isotope"]').click(function () {
            var isotope_array = $(this).attr("name").split("_");
            console.log(isotope_array);
            var adduct_dict = metabolite["adducts"][isotope_array[0]]["peaks"][isotope_array[1]];
            console.log(adduct_dict);
            $("#distribution_table tbody").html("");
            $("#dm_adduct").html(adduct_dict["type"]);
            for (i in adduct_dict["isotopic_distribution"]) {
                var spectra = adduct_dict["isotopic_distribution"][i];
                $("#distribution_table tbody").append("<tr><td>" + spectra[0].toFixed(4) + "</td><td>" + spectra[1].toFixed(2) + "</td></tr>")

            }
            $("#distribution_modal").modal("toggle");
        });
    });

}