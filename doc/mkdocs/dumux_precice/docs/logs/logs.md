# Overview over the last meetings

**Next meeting**: 2019-04-11

## 2019-04-03

Participants: Kilian, Ned, Alex

### Contents

- Discussion on the test case and how to set it up. 
- Some discussion about tools (`MkDocs`, `pandoc`, `sphinx` etc.). 
- Alex should come to one of the DuMuX meetings some time. Maybe present some of the tools.


### Outcome/Plans: 

- Create GitLab repository
- Set up [conjugate heat transfer test case](../dumux_test_case.md#conjugate-heat-transfer) in DuMuX. This will be done twice:
    1. Monolithic setup with DuMuX
    1. Two seperate DuMuX solvers that will be coupled with preCICE afterwards.
- Alex should push his notes about the `dumux-lecture` to the repository. 
- Checkpointing of the solution for preCICE has to be done "by hand". There is no functionality for this currently in DuMuX.
- Integrate the preCICE-DuMuX-adapter in the DuMuX repository. This keeps it close to potential users and it will not be forgotten in the preCICE repository.