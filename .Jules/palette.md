## 2024-04-18 - [SVG Data Legend]
**Learning:** Adding a static visual legend to generated SVGs (using reused styles and relative positioning) significantly improves comprehensibility for sighted users without affecting screen readers (which rely on existing `title` tags).
**Action:** Always include an explicit legend for data visualizations that use color or shape coding, even if `title` attributes exist.

## 2024-04-18 - [SVG Keyboard Accessibility]
**Learning:** Adding `tabindex="0"` and `:focus` styles to SVG elements (paths, shapes) makes them keyboard-accessible, allowing users to tab through data points and read titles/descriptions that would otherwise only be available on hover.
**Action:** Always ensure interactive or informative SVG elements are keyboard focusable and have visible focus indicators.

## 2024-04-18 - [Data-Rich Tooltips]
**Learning:** Embedding specific data values (like exact coordinates) into SVG element titles transforms them from generic labels into precise data inspection tools, benefiting both mouse users (hover) and screen reader users.
**Action:** Populate `title` attributes with dynamic, rounded data values rather than static text whenever possible.
