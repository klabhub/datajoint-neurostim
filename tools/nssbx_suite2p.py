"""Compatibility adapter for calling old and new Suite2p APIs from MATLAB."""

import suite2p


def _flatten_settings(settings, options):
    """Add nested Suite2p settings to the legacy flat options dictionary."""
    for key, value in settings.items():
        if isinstance(value, dict):
            _flatten_settings(value, options)
        else:
            options[key] = value


def default_options():
    """Return flat defaults for both the old and new Suite2p APIs."""
    if (hasattr(suite2p, "default_settings") and
            hasattr(suite2p, "default_db")):
        options = suite2p.default_db()
        _flatten_settings(suite2p.default_settings(), options)
        # Preserve the value type expected by the existing MATLAB code.
        if options.get("fast_disk") is None:
            options["fast_disk"] = []
        return options
    if hasattr(suite2p, "default_ops"):
        return suite2p.default_ops()
    raise AttributeError(
        "Suite2p exposes neither default_settings() nor default_ops()."
    )


def _apply_flat_settings(target, values):
    """Move matching flat values into a possibly nested settings dictionary."""
    for key, current_value in target.items():
        if isinstance(current_value, dict):
            _apply_flat_settings(current_value, values)
        elif key in values:
            target[key] = values.pop(key)


def run_suite2p(options, database):
    """Run Suite2p using the installed API generation."""
    if (hasattr(suite2p, "default_settings") and
            hasattr(suite2p, "default_db")):
        values = dict(options)
        values.update(database)

        new_database = suite2p.default_db()
        for key in tuple(new_database):
            if key in values:
                new_database[key] = values.pop(key)

        settings = suite2p.default_settings()
        _apply_flat_settings(settings, values)

        # Names used by older Suite2p releases.
        if "roidetect" in values:
            settings["run"]["do_detection"] = values.pop("roidetect")
        if "spikedetect" in values:
            settings["run"]["do_deconvolution"] = values.pop("spikedetect")

        if values:
            unknown = ", ".join(sorted(values))
            raise KeyError(f"Unknown Suite2p options: {unknown}")

        return suite2p.run_s2p(db=new_database, settings=settings)

    if hasattr(suite2p, "default_ops"):
        return suite2p.run_s2p(ops=options, db=database)

    raise AttributeError(
        "Suite2p exposes neither default_settings() nor default_ops()."
    )

