function array_to_text(a) {
        if (a.length != 0) {
            var text = "";
            for (i = 0; i < a.length; i++) {
                text += a[i] + "; ";
            }
        }
        else {
            var text = "None";
        }
        return text
}

function text_to_none(t) {
    if (t != null) {
        return t
    }
    else {
        return "None"
    }
}

function basic_information(metabolite) {
    document.title = metabolite["Name"] +" : DIMEdb - Direct Infusion MEtabolite Database";
    $("#name").html(metabolite["Name"]);

    function fill_info(metabolite) {
        $("#synonyms").html(array_to_text(metabolite["Synonyms"]));
        $("#iupacname").html(text_to_none(metabolite["IUPAC Name"]));
        $("#molform").html(metabolite["Molecular Formula"].replace(/([0-9]+)/g, '<sub>$1</sub>'));
        $("#smiles").html(text_to_none(metabolite["Identifiers"]["SMILES"]));
        $("#inchi").html(text_to_none(metabolite["Identifiers"]["InChI"]));

    }

    fill_info(metabolite);

}

function sidebar(metabolite) {

    function get_structure(id) {
        $.ajax({
            url: getBaseURL() + "gen_structure/" + id,
            success: function (result) {
                $("#structure").attr("src", "data:image/png;base64," + result);
            }
        });
}

    function fill_identifiers(identifiers) {
        var ignore = ["SMILES", "InChi"];
        for (var prop in identifiers) {
            if ($.inArray(prop, ignore) != 0) {
                $("#id_list").append("<li class='list-group-item'>"+prop+"</li>");
                if (identifiers[prop] != null) {

                }
            }
        }
    }

    get_structure(metabolite["_id"]);
    fill_identifiers(metabolite["Identifiers"])

}

function modals(metabolite) {
    function molecular_properties(properties) {
        $("#hbd").html(text_to_none(properties["Hydrogen Bond Donors"]));
        $("#hba").html(text_to_none(properties["Hydrogen Bond Acceptors"]));
        $("#rotbonds").html(text_to_none(properties["Rotatable Bonds"]));
        $("#formalcharge").html(text_to_none(properties["Formal Charge"]));
        $("#logp").html(properties["logP"].toFixed(4));
        $("#polarsurfacearea").html(text_to_none(properties["Polar Surface Area"].toFixed(4)));
        $("#fractionofsp3").html(properties["Fraction of SP3 Carbon"].toFixed(4))
        $("#mrvalues").html(properties["MR Values"].toFixed(4))
        $("#rings").html(text_to_none(properties["Rings"]));
        $("#aromaticrings").html(text_to_none(properties["Aromatic Rings"]));
        $("#heavyatoms").html(text_to_none(properties["Heavy Atoms"]));
        $("#secondaryamines").html(text_to_none(properties["Secondary Amines"]));
        $("#etheroxygens").html(text_to_none(properties["Ether Oxygens"]));
        $("#hydrogroups").html(text_to_none(properties["Hydroxy Groups"]));
        $("#carbacids").html(text_to_none(properties["Carboxylic Acids"]));
    }

    molecular_properties(metabolite["Properties"]);
}

function render_metabolite_view(metabolite_id) {
    var metabolite = get_metabolite(metabolite_id);

    basic_information(metabolite);
    sidebar(metabolite);
    modals(metabolite)
}