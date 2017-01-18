/*
    New search method.

    == Standard text search option, keep it as it is - it's fine.

    == PPI / Neutral

    New table to be created

    | Name | Adduct | Adduct m/z | Neutral m/z | Difference |


 */

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
        $("#met_name_modal").html(metabolite["name"]);
        $("#get_json").attr("href", base_url + "api/metabolite/?id__exact=" + metabolite["id"]);
        $("#molecular_formula").html(metabolite["molecular_formula"].replace(/([0-9]+)/g, '<sub>$1</sub>'));
        $("#molecular_formula_search").attr("href", base_url + "api/metabolites/?molecular_formula__exact=" + metabolite["molecular_formula"])
        $("#accurate_mass").html(metabolite["accurate_mass"]);
        $("#num_atoms").html(metabolite["num_atoms"]);

        for (p in metabolite["pathways"]) {
            var pw = metabolite["pathways"][p];
            $("#pathway_listgroup").append("<li class='list-group-item'>"+
                pw["name"]
                +"</li>");
        }

        for (source in metabolite["sources"]) {
            if (metabolite["sources"][source]) {
                var id = metabolite["sources"][source];
                if (source == "kegg_id") {
                    var kegg_url = "http://www.genome.jp/dbget-bin/www_bget?compound+"+String(id)
                    $("#data_sources").append("<a href='"+kegg_url+"'><button class='btn btn-success btn-sm btn-space'><i class='glyphicon glyphicon-link'></i> KEGG</button></a>");
                }

                if (source == "chebi_id") {
                    var chebi_url = "http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:"+String(id)
                    $("#data_sources").append("<a href='"+chebi_url+"'><button class='btn btn-danger btn-sm btn-space'><i class='glyphicon glyphicon-link'></i> CHEBI</button></a>");

                }

                if (source == "pubchem_id") {
                    var pubchem_url = "https://pubchem.ncbi.nlm.nih.gov/compound/"+String(id)
                    $("#data_sources").append("<a href='"+pubchem_url+"'><button class='btn btn-warning btn-sm btn-space'><i class='glyphicon glyphicon-link'></i> PubChem</button></a>");

                }
            }
        }

        // $("#data_sources").html("<a href='http://www.hmdb.ca/metabolites/"+metabolite["source"]+"' target='_blank'><button class='btn btn-primary'>"+metabolite["source"]+"</button></a>")
        $("#smiles").html(metabolite["smiles"]);

        for (indx in metabolite["synonyms"]) {
            $("#synonyms").append(metabolite["synonyms"][indx] + "; ");
        }

        for (indx in metabolite["origins"]) {
            $("#origins").append(metabolite["origins"][indx] + "; ");
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
                    $("#"+result+"_adduct").append("<li class='list-group-item'><b>"+adduct["type"]+":</b> "+adduct["accurate_mass"].toFixed(4)+"<button id='isotope' name='"+result+"_"+adduct_idx+"' class='btn btn-primary btn-sm pull-right'><i class='glyphicon glyphicon glyphicon-stats'></i> Isotope</button><div class='clearfix'></div></li>");
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
            type: 'bar',
            marker: {
                color: 'rgba(0, 0, 0, 1)'
            }
        }];

        var layout = {
            xaxis: {
                title: 'Mass-to-ion (m/z)',
                showgrid: false,
                range: [Math.min.apply(Math, x_plot)-2, Math.max.apply(Math, x_plot)+2]
            },
            yaxis: {
                title: 'Relative Intensity (%)'
            },
            margin: {
                l: 50,
                r: 50,
                b: 50,
                t: 50,
                pad: 2
            },
            bargap : 0.99
        };

        Plotly.newPlot("distribution_chart", isotopic_data, layout, {displayModeBar: false, barwidth :10});
        $("#container").show();

        $('[id="isotope"]').click(function () {
            var x_plot = []
            var y_plot = []
            var isotope_array = $(this).attr("name").split("_");
            var adduct_dict = metabolite["adducts"][isotope_array[0]]["peaks"][isotope_array[1]];
            $("#distribution_table tbody").html("");
            $("#dm_adduct").html(adduct_dict["type"]);
            for (i in adduct_dict["isotopic_distribution"]) {
                var spectra = adduct_dict["isotopic_distribution"][i];
                $("#distribution_table tbody").append("<tr><td>" + spectra[0].toFixed(4) + "</td><td>" + spectra[1].toFixed(2) + "</td></tr>")
                x_plot.push(spectra[0]);
                y_plot.push(spectra[1]);
            }

            var modal_data = [{
            x: x_plot,
            y: y_plot,
            type: 'bar',
            marker: {
                color: 'rgba(0, 0, 0, 1)'
            }
            }];

            Plotly.newPlot("distribution_chart_modal", modal_data, layout, {displayModeBar: false});

            $("#distribution_modal").modal("toggle");
        });
    });

}