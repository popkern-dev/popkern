//! Known failure classes (the ratchet): `cell/oracle/op` → the task that
//! fixes it. A failure not listed fails the run; a listed class that a
//! run covering its cell does not show must be removed (it is fixed).
//! Filled by the baseline inventory (PK-1236); every entry is a filed
//! task, and the list ends empty at the foundations' close-out.

pub const KNOWN: &[(&str, &str)] = &[
    ("SxS-O/policy/Union", "PK-1261"),
    ("SxS-O/policy/Common", "PK-1261"),
    ("SxS-O/policy/AminusB", "PK-1261"),
    ("SxS-O/policy/BminusA", "PK-1261"),
    ("YxY-O/policy/Common", "PK-1262"),
    ("YxY-O/policy/AminusB", "PK-1262"),
];
